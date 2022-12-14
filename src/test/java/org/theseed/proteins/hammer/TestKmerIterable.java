/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * Test the SequenceKmerIterable
 */
class TestKmerIterable {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TestKmerIterable.class);


    @Test
    void testKmerIterable() throws FileNotFoundException {
        // We read all the sequences and test them one at a time.
        File contigFile = new File("data", "contigs.fa");
        try (FastaInputStream contigStream = new FastaInputStream(contigFile)) {
            for (Sequence seq : contigStream) {
                log.info("Processing {} length {}.", seq.getLabel(), seq.length());
                this.processSequence(seq);
            }
        }
    }

    public static class Counter {
        int count;

        public Counter() {
            count = 0;
        }

        public void increment() {
            count++;
        }

        public int get() {
            return this.count;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + this.count;
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof Counter)) {
                return false;
            }
            Counter other = (Counter) obj;
            if (this.count != other.count) {
                return false;
            }
            return true;
        }
    }

    /**
     * @param seq	sequence to test
     */
    final static int K = 20;
    private void processSequence(Sequence seq) {
        String seqString = seq.getSequence();
        var originalCounts = new HashMap<String, Counter>(16000);
        int len = seqString.length();
        int limit = len - K;
        for (int i = 0; i <= limit; i++) {
            String kmer = seqString.substring(i, i + K);
            countKmer(originalCounts, kmer);
        }
        // The above is the ground truth.  Now use the iterator.
        SequenceKmerIterable kmers = new SequenceKmerIterable(seqString, K);
        var testCounts = new ConcurrentHashMap<String, Counter>(16000);
        for (String kmer : kmers)
            countKmer(testCounts, kmer);
        compareMaps(seq.getLabel(), "Iterable", originalCounts, testCounts);
        // Next, use a stream.
        testCounts.clear();
        kmers.stream().forEach(x -> countKmer(testCounts, x));
        compareMaps(seq.getLabel(), "Stream", originalCounts, testCounts);
        // Finally, a parallel stream.
        testCounts.clear();
        kmers.parallelStream().forEach(x -> countKmer(testCounts, x));
        compareMaps(seq.getLabel(), "ParallelStream", originalCounts, testCounts);
    }

    /**
     * Compare two maps and assert the counts are equal.
     *
     * @param label				label for the sequence being processed
     * @param method			method of iteration
     * @param originalCounts	ground truth
     * @param testCounts		counts to test
     */
    private void compareMaps(String label, String method, HashMap<String, Counter> originalCounts,
            ConcurrentHashMap<String, Counter> testCounts) {
        // Insure both maps are the same size.
        String testSeriesName = String.format("%s using %s", label, method);
        log.info("Analyzing {}.", testSeriesName);
        assertThat(testSeriesName, testCounts.size(), equalTo(originalCounts.size()));
        for (var testEntry : testCounts.entrySet()) {
            String kmer = testEntry.getKey();
            Counter testCount = testEntry.getValue();
            String testName = String.format("%s using %s, kmer = %s", label, method, kmer);
            assertThat(testName, originalCounts.containsKey(kmer), equalTo(true));
            assertThat(testName, testCount.get(), equalTo(originalCounts.get(kmer).get()));
        }
    }

    /**
     * Count a kmer in a map.
     *
     * @param countMap	counter map
     * @param kmer		kmer to count
     */
    private void countKmer(Map<String, Counter> countMap, String kmer) {
        assertThat(kmer, kmer.length(), equalTo(K));
        Counter c = countMap.computeIfAbsent(kmer, x -> new Counter());
        c.increment();
    }


}
