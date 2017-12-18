package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.tools.copynumber.caller.SimpleCopyRatioCaller;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.CopyRatioSegmentCollection;

import java.io.File;

/**
 * Call copy-ratio segments as amplified, deleted, or copy-number neutral.
 *
 * <h3>Input</h3>
 *
 * <li>
 *     Copy-ratio segments file.
 *     This is a .cr.seg file produced as one of the outputs of {@link ModelSegments}, which
 *     is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *     a row specifying the column headers contained in {@link CopyRatioSegmentCollection.CopyRatioSegmentTableColumn},
 *     and the corresponding entry rows.
 * </li>
 *
 * <h3>Output</h3>
 *
 * <li>
 *     Called copy-ratio segments file.
 *     This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *     a row specifying the column headers contained in {@link CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn},
 *     and the corresponding entry rows.
 * </li>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 *     gatk CallCopyRatioSegments \
 *          -I tumor.cr.seg \
 *          -O tumor.called.seg
 * </pre>
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Call copy-ratio segments as amplified, deleted, or copy-number neutral.",
        oneLineSummary = "Call copy-ratio segments as amplified, deleted, or copy-number neutral.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class CallCopyRatioSegments extends CommandLineProgram {
    public static final String NEUTRAL_SEGMENT_COPY_RATIO_THRESHOLD_LONG_NAME = "neutral-segment-copy-ratio-threshold";
    public static final String OUTLIER_NEUTRAL_SEGMENT_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME = "outlier-neutral-segment-copy-ratio-z-score-threshold";
    public static final String CALLING_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME = "calling-copy-ratio-z-score-threshold";

    @Argument(
            doc = "Input file containing copy-ratio segments (.cr.seg output of ModelSegments).",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME
    )
    private File segmentsFile;

    @Argument(
            doc = "Output file for called segments.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outFile;

    @Argument(
            doc = "Threshold on non-log2 copy ratio used for determining copy-neutral segments.  " +
                    "If non-log2 copy ratio is within 1 +/- this threshold, a segment is considered copy-neutral.",
            fullName = NEUTRAL_SEGMENT_COPY_RATIO_THRESHOLD_LONG_NAME,
            optional = true
    )
    private double neutralSegmentCopyRatioThreshold = 0.1;

    @Argument(
            doc = "Threshold on z-score of non-log2 copy ratio used for determining outlier copy-neutral segments.  " +
                    "If non-log2 copy ratio z-score is above this threshold for a copy-neutral segment, " +
                    "it is considered an outlier and not used in the calculation of the length-weighted mean and standard deviation " +
                    "used for calling.",
            fullName = OUTLIER_NEUTRAL_SEGMENT_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double outlierNeutralSegmentCopyRatioZScoreThreshold = 2.;

    @Argument(
            doc = "Threshold on z-score of non-log2 copy ratio used for calling segments.",
            fullName = CALLING_COPY_RATIO_Z_SCORE_THRESHOLD_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double callingCopyRatioZScoreThreshold = 2.;

    @Override
    protected Object doWork() {
        final CopyRatioSegmentCollection copyRatioSegments = new CopyRatioSegmentCollection(segmentsFile);
        final CalledCopyRatioSegmentCollection calledCopyRatioSegments =
                new SimpleCopyRatioCaller(copyRatioSegments,
                        neutralSegmentCopyRatioThreshold, outlierNeutralSegmentCopyRatioZScoreThreshold, callingCopyRatioZScoreThreshold)
                        .makeCalls();
        calledCopyRatioSegments.write(outFile);

        return "SUCCESS";
    }
}
