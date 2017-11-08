package org.broadinstitute.hellbender.utils.runtime;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.IOException;
import java.io.OutputStream;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.Consumer;

/**
 * A service that can be used to write to a stream using a thread background thread and an executor service. This
 * is typically used to write items to a buffered stream that might block until the stream is consumed by a reader.
 * @param <T> Type of items to be written.
 */
public class AsynchronousStreamWriterService<T> {

    private static final Logger logger = LogManager.getLogger(AsynchronousStreamWriterService.class);

    final ExecutorService executorService;
    final OutputStream streamWriter;
    Future<Integer> previousBatch;

    /**
     * @param executorService executor service to be used to dispatch background tasks
     * @param streamWriter target stream to which items should be written
     */
    public AsynchronousStreamWriterService(
            final ExecutorService executorService,
            final OutputStream streamWriter)
    {
        Utils.nonNull(executorService);
        Utils.nonNull(streamWriter);

        this.streamWriter = streamWriter;
        this.executorService = executorService;
        previousBatch = null;
    }

    /**
     * Request that a batch of items be written to the stream on a background thread. Any previously requested batch
     * must have already been completed and retrieved via {@link #waitForPreviousBatchCompletion}.
     *
     * @param batchList a list of items to be written
     * @param batchSize number of items to write
     * @return {@code Future} representing this batch
     */
    public Future<Integer> startAsynchronousBatchWrite(final List<T> batchList, int batchSize) {
        Utils.nonNull(batchList);
        Utils.nonEmpty(batchList);

        if (previousBatch != null) {
            throw new IllegalStateException("Previous batch not yet complete");
        }

        previousBatch = executorService.submit(() -> {
            try {
                for (int i = 0; i < batchSize; i++) {
                    T element = batchList.get(i);
                    streamWriter.write(element.toString().getBytes());
                }
                // this can block, waiting for the stream to be consumed if its buffered
                streamWriter.flush();
                return batchSize; // return the number of items this batch was asked to write
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        return previousBatch;
    }

    /**
     * Waits for a batch that was previously initiated via {@link #startAsynchronousBatchWrite(List, int)}}
     * to complete, flushes the target stream and returns the corresponding completed Future. The Future representing
     * a given batch can only be obtained via this method once. If no work is outstanding, and/or the previous batch
     * has already been retrieved, null is returned.
     * @param timeout maximum time to wait for a response
     * @param timeoutUnit units of {@code timeout}
     * @return returns null if no previous work to complete, otherwise a completed Future
     */
    public Future<Integer> waitForPreviousBatchCompletion(long timeout, TimeUnit timeoutUnit) {
        final Future<Integer> lastCompleteBatch = previousBatch;
        if (previousBatch != null) {
            try {
                try {
                    previousBatch.get(timeout, timeoutUnit);
                } catch (ExecutionException | InterruptedException e) {
                    throw new GATKException("Interrupted during background stream write");
                } catch (TimeoutException e) {
                    throw new GATKException("Timeout waiting for background stream write to complete");
                }
                streamWriter.flush();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            previousBatch = null;
        }
        return lastCompleteBatch;
    }

    /**
     * Terminate the async writer, cancelling any outstanding work.
     * @return true
     */
    public boolean terminate() {
        boolean isCancelled = true;
        if (previousBatch != null) {
            logger.warn("Cancelling outstanding asynchronous writing");
            isCancelled = previousBatch.cancel(true);
        }
        previousBatch = null;
        return isCancelled;
    }

    protected void finalize() throws Throwable {
        terminate();
    }

}
