/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * This sequence manager reads all the sequences at the beginning and loads them into memory.  They
 * are then served by list iterators.
 *
 * @author Bruce Parrello
 *
 */
public class CacheSequenceManager extends SequenceManager {

    // FIELDS
    /** list of sequences read */
    private List<Sequence> seqs;

    /**
     * This is the class to create an iterator.  Since there is no open file, it's pretty simple.
     */
    protected class Iter extends SequenceManager.ResourceIter {

        /** iterator through the list */
        private Iterator<Sequence> iter;

        /**
         * Construct the iterator.  We simply pass through an iterator on the sequence list.
         */
        public Iter() {
            this.iter = CacheSequenceManager.this.seqs.iterator();
        }

        @Override
        public boolean hasNext() {
            return this.iter.hasNext();
        }

        @Override
        public Sequence next() {
            return this.iter.next();
        }

        @Override
        public void close() {
        }



    }

    /**
     * Construct a new cached-sequence manager.
     *
     * @param file	name of a FASTA file containing the sequences
     * @throws IOException
     */
    public CacheSequenceManager(File file) throws IOException {
        super(file);
        // Create the list.  We estimate the list size using mean contig lengths, which will be grossly low for
        // feature sequences, but will still reduce the number of attempts to grow the list.
        int estimatedSize = (int) (file.length() / 8000);
        this.seqs = new ArrayList<Sequence>(estimatedSize);
        log.info("Reading FASTA file {} into cache.", file);
        long len = 0;
        int count = 0;
        try (var inStream = new FastaInputStream(file)) {
            for (Sequence seq : inStream) {
                this.seqs.add(seq);
                this.seqs.add(seq.reverse());
                len += seq.length();
                count++;
                if (log.isInfoEnabled() && count % 5000 == 0)
                    log.info("{} sequences read from file.", count);
            }
        }
        log.info("{} sequences with {} letters read from file.", count, len);
    }

    @Override
    public ResourceIter newIterator() {
        // Return an iterator through the list.
        return this.new Iter();
    }

}
