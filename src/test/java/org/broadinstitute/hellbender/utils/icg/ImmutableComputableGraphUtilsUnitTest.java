package org.broadinstitute.hellbender.utils.icg;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.annotations.Test;

/**
 * Unit tests for {@link ImmutableComputableGraphUtils}.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public class ImmutableComputableGraphUtilsUnitTest extends GATKBaseTest {

    private static final CacheNode.NodeKey X_KEY = new CacheNode.NodeKey("x");

    @Test(expectedExceptions = ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder.DuplicateNodeKeyException.class)
    public void testDuplicatePrimitiveNode() {
        ImmutableComputableGraph.builder()
                .primitiveNode(X_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .primitiveNode(X_KEY, new CacheNode.NodeTag[] {}, new DuplicableNDArray())
                .build();
    }

    @Test(expectedExceptions = ImmutableComputableGraphUtils.ImmutableComputableGraphBuilder.DuplicateNodeKeyException.class)
    public void testDuplicateComputableNode() {
        ImmutableComputableGraph.builder()
                .computableNode(X_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {}, null, true)
                .computableNode(X_KEY, new CacheNode.NodeTag[] {}, new CacheNode.NodeKey[] {}, null, true)
                .build();
    }
}
