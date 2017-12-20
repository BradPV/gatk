package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.tools.walkers.validation.basicshortmutpileup.BetaBinomialDistribution;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;

import java.io.File;
import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Given pileup data from {@link GetPileupSummaries}, calculates the fraction of reads coming from cross-sample contamination.
 *
 * <p>
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 *
 * <p>This tool and GetPileupSummaries together replace GATK3's ContEst.</p>
 *
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 * <h3>Example</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Calculate the fraction of reads coming from cross-sample contamination",
        oneLineSummary = "Calculate the fraction of reads coming from cross-sample contamination",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    private static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    private static final int MAX_CHANGEPOINTS_PER_GENOME = 100;
    private static final int MAX_CHANGEPOINTS_PER_CHROMOSOME = 10;
    private static final int MIN_SITES_PER_SEGMENT = 5;

    // our analysis only cares about hom alt and het sites, so we throw away hom refs with a very conservative heuristic
    private static final double ALT_FRACTION_OF_DEFINITE_HOM_REF = 0.05;

    private static final double KERNEL_SEGMENTER_LINEAR_COST = 1.0;
    private static final double KERNEL_SEGMENTER_LOG_LINEAR_COST = 1.0;
    private static final int KERNEL_SEGMENTER_DIMENSION = 100;
    private static final int POINTS_PER_SEGMENTATION_WINDOW = 50;

    private static final double DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD = 1.0/2;
    private static final double DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD = 2.0;

    // shape parameters for beta binomial distributions of hom alt and het read counts
    // the het distribution is peaked at an alt fraction of 0.5 and falls to 1/10 of its max around 0.35 and 0.65
    private static final double HET_ALPHA = 30.0;
    private static final double HET_BETA = 30.0;
    private double HOM_ALT_ALPHA = 120.0;
    private static final double HOM_ALT_BETA = 1.0;

    private static final List<Double> HET_ALPHA_BETA_TO_TRY = Arrays.asList(5., 10., 25., 50., 100., 200., 500.);
    private static final List<Double> HOM_ALT_ALPHA_TO_TRY = Arrays.asList(5., 10., 25., 50., 100., 200., 500.);


    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table")
    private File inputPileupSummariesTable;

    @Argument(fullName = "matchedNormal",
            shortName = "matched",
            doc="The matched normal input table", optional = true)
    private File matchedPileupSummariesTable = null;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table")
    private final File outputTable = null;

    @Argument(fullName= "lowCoverageRatioThreshold",
            doc="The minimum coverage relative to the median.", optional = true)
    private final double lowCoverageRatioThreshold = DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD;

    @Argument(fullName= "highCoverageRatioThreshold",
            doc="The maximum coverage relative to the median.", optional = true)
    private final double highCoverageRatioThreshold = DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD;

    private static final double SEGMENTATION_KERNEL_VARIANCE = 0.025;
    
    private static final BiFunction<PileupSummary, PileupSummary, Double> SEGMENTATION_KERNEL = (ps1, ps2) -> {
        final double maf1 = FastMath.min(ps1.getAltFraction(), 1 - ps1.getAltFraction());
        final double maf2 = FastMath.min(ps2.getAltFraction(), 1 - ps2.getAltFraction());
        return FastMath.exp(-MathUtils.square(maf1 - maf2)/(2 * SEGMENTATION_KERNEL_VARIANCE));
    };

    @Override
    public Object doWork() {
        final List<PileupSummary> pileupSummaries = filterSites(PileupSummary.readPileupSummaries(inputPileupSummariesTable));
        final List<PileupSummary> matchedPileupSummaries = matchedPileupSummariesTable == null ? null
                : filterSites(PileupSummary.readPileupSummaries(matchedPileupSummariesTable));

        // used the matched normal to genotype if available
        final List<PileupSummary> pileupSummariesForGenotyping = matchedPileupSummaries == null ? pileupSummaries : matchedPileupSummaries;
        final double expectedHomAltCount = pileupSummariesForGenotyping.stream()
                .mapToDouble(PileupSummary::getAlleleFrequency)
                .map(MathUtils::square)
                .sum();

        int iteration = 0;
        double currentContamination = 0.0;
        while (true) {
            final double oldContamination = currentContamination;
            final List<PileupSummary> homAltSites = findHomAltSites(pileupSummariesForGenotyping);
            final List<PileupSummary> homAltPileups = matchedPileupSummaries == null ? homAltSites : getCorrespondingHomAltSites(homAltSites, pileupSummaries);
            final Pair<Double, Double> contaminationAndError = calculateContamination(homAltPileups, true);
            currentContamination = contaminationAndError.getLeft();

            if (homAltSites.size() > expectedHomAltCount) {
                HOM_ALT_ALPHA *= 1.3;
            } else {
                HOM_ALT_ALPHA /= 1.3;
            }


            if (++iteration > 10 || Math.abs(currentContamination - oldContamination) < 0.001) {
                final double error = contaminationAndError.getRight();
                ContaminationRecord.writeContaminationTable(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), currentContamination, error)), outputTable);
                break;
            }
        }

        return "SUCCESS";
    }

    private static List<PileupSummary> getCorrespondingHomAltSites(final List<PileupSummary> homAltsInMatchedNormal, final List<PileupSummary> allSites) {
        final OverlapDetector<PileupSummary> homAltsInMatchedNormalOverlapDetector = OverlapDetector.create(homAltsInMatchedNormal);
        return allSites.stream().filter(homAltsInMatchedNormalOverlapDetector::overlapsAny).collect(Collectors.toList());
    }

    private List<PileupSummary> findHomAltSites(List<PileupSummary> pileupSummaries) {
        final List<IndexRange> segments = findAllelicCopyNumberSegments(pileupSummaries)
                .stream().filter(s -> s.size() >= MIN_SITES_PER_SEGMENT).collect(Collectors.toList());

        // test
        final List<List<PileupSummary>> segmentPileups = segments.stream()
                .map(segment -> IntStream.range(segment.from, segment.to).mapToObj(pileupSummaries::get).collect(Collectors.toList()))
                .collect(Collectors.toList());
        final double x = genomeLogLikelihood(segmentPileups);
        // test

        final List<List<PileupSummary>> homAltSitesBySegment = findHomAltSitesBySegment(pileupSummaries, segments);

        final double[] contaminationEstimatesBySegment = homAltSitesBySegment.stream()
                .mapToDouble(segmentSites -> calculateContamination(segmentSites, false).getLeft())
                .toArray();
        final double[] nonZeroContaminationEstimatesBySegment = Arrays.stream(contaminationEstimatesBySegment)
                .filter(c -> c > 0).toArray();
        final double medianSegmentEstimate = new Median().evaluate(nonZeroContaminationEstimatesBySegment);

        // filter out segments with very high contamination, which is likely due to loss of heterozygosity
        final List<List<PileupSummary>> believableHomAltSitesBySegment = IntStream.range(0, homAltSitesBySegment.size())
                .filter(n -> contaminationEstimatesBySegment[n] < 2 * medianSegmentEstimate || contaminationEstimatesBySegment[n] < medianSegmentEstimate + 0.05)
                .mapToObj(homAltSitesBySegment::get)
                .collect(Collectors.toList());

        return believableHomAltSitesBySegment.stream().flatMap(List::stream).collect(Collectors.toList());
    }


    private List<List<PileupSummary>> findHomAltSitesBySegment(List<PileupSummary> pileupSummaries, List<IndexRange> segments) {
        final List<List<PileupSummary>> homAltSitesBySegment = new ArrayList<>();
        for (final IndexRange segment : segments) {
            logger.info(String.format("Considering segment with %d sites from %s:%d to %s:%d", segment.size(),
                    pileupSummaries.get(segment.from).getContig(), pileupSummaries.get(segment.from).getStart(),
                    pileupSummaries.get(segment.to - 1).getContig(), pileupSummaries.get(segment.to - 1).getStart()));

            final List<PileupSummary> segmentSites = IntStream.range(segment.from, segment.to).mapToObj(pileupSummaries::get).collect(Collectors.toList());
            final List<PileupSummary> homAltSitesInSegment = findHomAltSitesInSegment(segmentSites);

            homAltSitesBySegment.add(homAltSitesInSegment);
        }
        return homAltSitesBySegment;
    }

    private List<IndexRange> findAllelicCopyNumberSegments(List<PileupSummary> pileupSummaries) {
        final List<Integer> changepoints = new ArrayList<>();

        // when the kernel segmenter finds a changepoint at index n, that means index n belongs to the left segment, which goes
        // against the usual end-exclusive intervals of IndexRange etc.  This explains adding in the first changepoint of -1
        // instead of 0 and all the "changepoint + 1" constructions below
        changepoints.add(-1);
        changepoints.addAll(getChangepoints(pileupSummaries));
        changepoints.add(pileupSummaries.size()-1);
        return IntStream.range(0, changepoints.size() - 1)
                .mapToObj(n -> new IndexRange(changepoints.get(n) + 1, changepoints.get(n+1) + 1))
                .collect(Collectors.toList());
    }

    private static Pair<Double, Double> calculateContamination(List<PileupSummary> homAltSites, final boolean doLogging) {
        if (homAltSites.isEmpty()) {
            logger.warn("No hom alt sites found!  Perhaps GetPileupSummaries was run on too small of an interval, or perhaps the sample was extremely inbred or haploid.");
            return Pair.of(0.0, 1.0);
        }

        final long totalReadCount = homAltSites.stream().mapToLong(PileupSummary::getTotalCount).sum();
        final long totalRefCount = homAltSites.stream().mapToLong(PileupSummary::getRefCount).sum();

        // if eg ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final long errorRefCount = homAltSites.stream().mapToLong(PileupSummary::getOtherAltCount).sum() / 2;
        final long contaminationRefCount = Math.max(totalRefCount - errorRefCount, 0);
        final double totalDepthWeightedByRefFrequency = homAltSites.stream()
                .mapToDouble(ps -> ps.getTotalCount() * (1 - ps.getAlleleFrequency()))
                .sum();
        final double contamination = contaminationRefCount / totalDepthWeightedByRefFrequency;
        final double standardError = Math.sqrt(contamination / totalDepthWeightedByRefFrequency);

        if (doLogging) {
            logger.info(String.format("In %d homozygous variant sites we find %d reference reads due to contamination and %d" +
                    " due to to sequencing error out of a total %d reads.", homAltSites.size(), contaminationRefCount, errorRefCount, totalReadCount));
            logger.info(String.format("Based on population data, we would expect %d reference reads in a contaminant with equal depths at these sites.", (long) totalDepthWeightedByRefFrequency));
            logger.info(String.format("Therefore, we estimate a contamination of %.3f.", contamination));
            logger.info(String.format("The error bars on this estimate are %.5f.", standardError));
        }
        return Pair.of(contamination, standardError);
    }

    private List<PileupSummary> filterSites(final List<PileupSummary> allSites) {
        final double[] coverage = allSites.stream().mapToDouble(PileupSummary::getTotalCount).toArray();
        final double medianCoverage = new Median().evaluate(coverage);
        final double lowCoverageThreshold = medianCoverage * lowCoverageRatioThreshold;
        final double highCoverageThreshold = medianCoverage * highCoverageRatioThreshold;
        return allSites.stream()
                .filter(ps -> ps.getTotalCount() > lowCoverageThreshold && ps.getTotalCount() < highCoverageThreshold)
                .filter(ps -> ps.getAltFraction() > ALT_FRACTION_OF_DEFINITE_HOM_REF)
                .collect(Collectors.toList());
    }

    private static List<Integer> getChangepoints(final List<PileupSummary> pileups) {
        final int numChromosomes = (int) pileups.stream().map(PileupSummary::getContig).distinct().count();
        final int numChangepoints = FastMath.min(numChromosomes * MAX_CHANGEPOINTS_PER_CHROMOSOME, MAX_CHANGEPOINTS_PER_GENOME);
        final KernelSegmenter<PileupSummary> segmenter = new KernelSegmenter<>(pileups);
        return segmenter.findChangepoints(numChangepoints, SEGMENTATION_KERNEL, KERNEL_SEGMENTER_DIMENSION,
                Arrays.asList(POINTS_PER_SEGMENTATION_WINDOW), KERNEL_SEGMENTER_LINEAR_COST, KERNEL_SEGMENTER_LOG_LINEAR_COST, KernelSegmenter.ChangepointSortOrder.INDEX);
    }


    private List<PileupSummary> findHomAltSitesInSegment(final List<PileupSummary> sites) {
        return sites.stream().filter(site -> homAltProbability(site) > 0.5).collect(Collectors.toList());
    }

    private double homAltProbability(final PileupSummary site) {
        final double alleleFrequency = site.getAlleleFrequency();
        final double homAltPrior = MathUtils.square(alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homRefPrior = MathUtils.square(1 - alleleFrequency);

        final int altCount = site.getAltCount();
        final int totalCount = altCount + site.getRefCount();

        final double homAltLikelihood = new BetaBinomialDistribution(null, HOM_ALT_ALPHA, HOM_ALT_BETA, totalCount).probability(altCount);
        final double hetLikelihood = new BetaBinomialDistribution(null, HET_ALPHA, HET_BETA, totalCount).probability(altCount);

        final double unnormalizedHomAltProbability = homAltPrior * homAltLikelihood;
        final double unnormalizedHetProbability = hetPrior * hetLikelihood;

        final double result = unnormalizedHomAltProbability / (unnormalizedHetProbability + unnormalizedHomAltProbability);

        return result;

    }

    /**
     * Given the beta binomial distributions with shapes
     * het: alpha = beta = hetaAlphaBeta
     * homAlt: alpha = homAltAlpha, beta = 1
     * homRef: alpha = 1, beta = homAltAlpha
     * find the log likelihood of a site
     */
    private double siteLogLikelihood(final PileupSummary site, final double hetAlphaBeta, final double homAltAlpha) {
        final double alleleFrequency = site.getAlleleFrequency();
        final double homAltPrior = MathUtils.square(alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);
        final double homRefPrior = 0.1; //MathUtils.square(1 - alleleFrequency);

        final int altCount = site.getAltCount();
        final int totalCount = altCount + site.getRefCount();

        // debug
        if (altCount < totalCount / 3) {
            return 0;
        }
        // debug

        final double homAltLikelihood = new BetaBinomialDistribution(null, homAltAlpha, 1, totalCount).probability(altCount);
        final double homRefLikelihood = new BetaBinomialDistribution(null, 1, homAltAlpha, totalCount).probability(altCount);
        final double hetLikelihood = new BetaBinomialDistribution(null, hetAlphaBeta, hetAlphaBeta, totalCount).probability(altCount);

        return FastMath.log(homAltPrior * homAltLikelihood + hetPrior * hetLikelihood + homRefPrior * homRefLikelihood);
    }

    /**
     * the log likelihood of all sites in a segment
     */
    private double segmentLogLikelihood(final List<PileupSummary> sites, final double hetAlphaBeta, final double homAltAlpha) {
        return sites.stream().mapToDouble(site -> siteLogLikelihood(site, hetAlphaBeta, homAltAlpha)).sum();
    }

    private double segmentLogLikelihood(final List<PileupSummary> sites, final double homAltAlpha) {
        // initialize best index as the one corresponding to too much variance in alt fraction i.e. loss of heterozygosity
        int bestIndex = 0;
        double bestLogLikelihood = Double.NEGATIVE_INFINITY;
        for (int n = 0; n < HET_ALPHA_BETA_TO_TRY.size(); n++) {
            final double hetAlphaBeta = HET_ALPHA_BETA_TO_TRY.get(n);
            final double logLikelihood = segmentLogLikelihood(sites, hetAlphaBeta, homAltAlpha);
            if (logLikelihood > bestLogLikelihood) {
                bestLogLikelihood = logLikelihood;
                bestIndex = n;
            }
        }
        return bestLogLikelihood;
    }

    private double genomeLogLikelihood(final List<List<PileupSummary>> segments, final double homAltAlpha) {
        return segments.stream().mapToDouble(segment -> segmentLogLikelihood(segment, homAltAlpha)).sum();
    }

    private double genomeLogLikelihood(final List<List<PileupSummary>> segments) {
        int bestIndex = 0;
        double bestLogLikelihood = Double.NEGATIVE_INFINITY;
        for (int n = 0; n < HOM_ALT_ALPHA_TO_TRY.size(); n++) {
            final double homAltAlpha = HOM_ALT_ALPHA_TO_TRY.get(n);
            final double logLikelihood = genomeLogLikelihood(segments, homAltAlpha);
            if (logLikelihood > bestLogLikelihood) {
                bestLogLikelihood = logLikelihood;
                bestIndex = n;
            }
        }
        return bestLogLikelihood;
    }
}
