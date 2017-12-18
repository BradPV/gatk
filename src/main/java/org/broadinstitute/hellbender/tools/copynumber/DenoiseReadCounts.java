package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hdf5.HDF5Library;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.denoising.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.utils.IntervalMergingRule;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;

/**
 * Denoise read counts to produce a denoised copy-ratio profile.
 *
 * <p>
 *     Typically, a panel of normals produced by {@link CreateReadCountPanelOfNormals} is provided as input.
 *     The input counts are then standardized by 1) transforming to fractional coverage,
 *     2) performing optional explicit GC-bias correction (if the panel contains GC-content annotated intervals),
 *     3) filtering intervals to those contained in the panel, 4) dividing by interval medians contained in the panel,
 *     5) dividing by the sample median, and 6) transforming to log-2 copy ratio.  The result is then denoised by
 *     subtracting the projection onto the specified number of principal components from the panel.
 * </p>
 *
 * <p>
 *     If no panel is provided, then the input counts are instead standardized by 1) transforming to fractional coverage,
 *     2) performing optional explicit GC-bias correction (if GC-content annotated intervals are provided),
 *     3) dividing by the sample median, and 4) transforming to log-2 copy ratio.  The denoised result is taken to be
 *     identical to the standardized result.
 * </p>
 *
 * <h3>Input</h3>
 *
 * <li>
 *     Counts file.
 *     This is the TSV or HDF5 output of {@link CollectFragmentCounts}.
 * </li>
 * <li>
 *     (Optional) Panel-of-normals file.
 *     This is the output of {@link CreateReadCountPanelOfNormals}.  If provided,
 *     it will be used to standardize and denoise the input counts.  This may include explicit GC-bias correction
 *     if annotated intervals were used to create the panel.
 * </li>
 * <li>
 *     (Optional) GC-content annotated-intervals file.
 *     This can be provided in place of a panel of normals to perform explicit GC-bias correction.
 * </li>
 *
 * <h3>Output</h3>
 *
 * <li>
 *     Standardized copy-ratio profile file.
 *     This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *     a row specifying the column headers contained in {@link CopyRatioCollection.CopyRatioTableColumn},
 *     and the corresponding entry rows.
 * </li>
 * <li>
 *     Denoised copy-ratio profile file.
 *     This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *     a row specifying the column headers contained in {@link CopyRatioCollection.CopyRatioTableColumn},
 *     and the corresponding entry rows.
 * </li>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --readCountPanelOfNormals panel_of_normals.pon.hdf5 \
 *          --standardizedCopyRatios sample.standardizedCR.tsv \
 *          --denoisedCopyRatios sample.denoisedCR.tsv
 * </pre>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --annotated-intervals annotated_intervals.tsv \
 *          --standardizedCopyRatios sample.standardizedCR.tsv \
 *          --denoisedCopyRatios sample.denoisedCR.tsv
 * </pre>
 *
 * <pre>
 *     gatk DenoiseReadCounts \
 *          -I sample.counts.hdf5 \
 *          --standardizedCopyRatios sample.standardizedCR.tsv \
 *          --denoisedCopyRatios sample.denoisedCR.tsv
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoise read counts using a panel of normals.",
        oneLineSummary = "Denoise read counts using a panel of normals.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class DenoiseReadCounts extends CommandLineProgram {
    @Argument(
            doc = "Input TSV or HDF5 read-count file containing integer read counts in genomic intervals for a single case sample.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private File inputReadCountFile;

    @Argument(
            doc = "Input HDF5 file containing the panel of normals (output of CreateReadCountPanelOfNormals).",
            fullName = CopyNumberStandardArgument.COUNT_PANEL_OF_NORMALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputPanelOfNormalsFile = null;

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals (output of AnnotateIntervals).  " +
                    "Intervals must be identical to and in the same order as those in the input read-count file.  " +
                    "If a panel of normals is provided, this input will be ignored.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File annotatedIntervalsFile = null;

    @Argument(
            doc = "Output file for standardized copy-ratio profile.  GC-bias correction will be performed if annotations for GC content are provided.",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File standardizedCopyRatiosFile;

    @Argument(
            doc = "Output file for denoised copy-ratio profile.",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME
    )
    private File denoisedCopyRatiosFile;

    @Argument(
            doc = "Number of eigensamples to use for denoising.  " +
                    "If not specified or if the number of eigensamples available in the panel of normals " +
                    "is smaller than this, all eigensamples will be used.",
            fullName = CopyNumberStandardArgument.NUMBER_OF_EIGENSAMPLES_LONG_NAME,
            optional = true
    )
    private Integer numEigensamplesRequested = null;

    @Override
    protected Object doWork() {
        if (!new HDF5Library().load(null)) { //Note: passing null means using the default temp dir.
            throw new UserException.HardwareFeatureException("Cannot load the required HDF5 library. " +
                    "HDF5 is currently supported on x86-64 architecture and Linux or OSX systems.");
        }
        Utils.validateArg(numEigensamplesRequested == null || numEigensamplesRequested > 0,
                "Number of eigensamples to use for denoising must be non-negative.");

        IOUtils.canReadFile(inputReadCountFile);
        logger.info(String.format("Reading read-count file (%s)...", inputReadCountFile));
        final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);

        if (inputPanelOfNormalsFile != null) {  //denoise using panel of normals
            IOUtils.canReadFile(inputPanelOfNormalsFile);
            try (final HDF5File hdf5PanelOfNormalsFile = new HDF5File(inputPanelOfNormalsFile)) {  //HDF5File implements AutoCloseable
                final SVDReadCountPanelOfNormals panelOfNormals = HDF5SVDReadCountPanelOfNormals.read(hdf5PanelOfNormalsFile);

                if (annotatedIntervalsFile != null) {
                    logger.warn("Panel of normals was provided; ignoring input GC-content annotations...");
                }

                //perform denoising and write result
                final int numEigensamples =
                        numEigensamplesRequested == null ?
                                panelOfNormals.getNumEigensamples() :
                                Math.min(panelOfNormals.getNumEigensamples(), this.numEigensamplesRequested);
                if (numEigensamplesRequested != null && numEigensamples < numEigensamplesRequested) {
                    logger.warn(String.format("%d eigensamples were requested but only %d are available in the panel of normals...",
                            numEigensamplesRequested, numEigensamples));
                }
                final SVDDenoisedCopyRatioResult denoisedCopyRatioResult = panelOfNormals.denoise(readCounts, numEigensamples);

                logger.info("Writing standardized and denoised copy-ratio profiles...");
                denoisedCopyRatioResult.write(standardizedCopyRatiosFile, denoisedCopyRatiosFile);
            }
        } else {    //standardize and perform optional GC-bias correction
            //get GC content (null if not provided)
            final double[] intervalGCContent = GCBiasCorrector.validateIntervalGCContent(
                    readCounts.getMetadata().getSequenceDictionary(), readCounts.getIntervals(), annotatedIntervalsFile);

            if (intervalGCContent == null) {
                logger.warn("Neither a panel of normals nor GC-content annotations were provided, so only standardization will be performed...");
            }

            final RealMatrix standardizedCopyRatioValues = SVDDenoisingUtils.preprocessAndStandardizeSample(readCounts.getCounts(), intervalGCContent);

            //construct a result with denoised profile identical to standardized profile
            final SVDDenoisedCopyRatioResult standardizedResult = new SVDDenoisedCopyRatioResult(
                    readCounts.getMetadata(),
                    readCounts.getIntervals(),
                    standardizedCopyRatioValues,
                    standardizedCopyRatioValues);
            standardizedResult.write(standardizedCopyRatiosFile, denoisedCopyRatiosFile);
        }

        logger.info("Read counts successfully denoised.");

        return "SUCCESS";
    }
}
