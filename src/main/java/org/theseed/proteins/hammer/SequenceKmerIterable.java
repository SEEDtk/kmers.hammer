/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Iterator;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.theseed.sequence.Sequence;

/**
 * This object presents the kmers in a sequence as an iterable.  It also supports streams and parallel streams.  Unlike
 * KmerSeries, it has no support for a set of sequences:  only a single sequence is allowed.  This enables us to
 * do the splitting more efficiently.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceKmerIterable implements Iterable<String> {

    // FIELDS
    /** sequence string over which to iterate */
    private final String seqString;
    /** kmer size */
    private final int kmerSize;

    /**
     * This is an iterator for the kmers in the sequence.
     */
    protected class Iter implements Iterator<String> {

        /** current position in the sequence */
        private int pos;
        /** position past last permissible point in the sequence */
        private int limit;

        /**
         * Construct the iterator.
         */
        public Iter() {
            this.pos = 0;
            this.limit = SequenceKmerIterable.this.getLimit();
        }

        @Override
        public boolean hasNext() {
            return (pos < limit);
        }

        @Override
        public String next() {
            String retVal = SequenceKmerIterable.this.getKmer(pos);
            pos++;
            return retVal;
        }

    }

    /**
     * This is a spliterator for all the kmers in the sequence.  The spliterator provides support
     * for dividing up the iteration task into sub-tasks for parallelism purposes.  As such, the
     * spliterator covers a subsequence, and the first one is created covering the whole sequence.
     * A returned kmer may extend past the end of the chosen subsequence, but it can never extend
     * past the end of the whole sequence.  The first subsequence must place its initial limit so
     * the latter problem does not occur.
     */
    protected class Splitter implements Spliterator<String> {

        /** initial position in the sequence */
        private int start;
        /** current position in the sequence */
        private int pos;
        /** position past last permissible point in the sequence */
        private int limit;

        /**
         * Create a spliterator over the specified portion of the sequence.
         *
         * @param i1	position of first kmer
         * @param i2	position past start of last kmer
         */
        public Splitter(int i1, int i2) {
            this.pos = i1;
            this.start = i1;
            this.limit = i2;
        }

        @Override
        public boolean tryAdvance(Consumer<? super String> action) {
            boolean retVal;
            if (pos >= limit)
                retVal = false;
            else {
                action.accept(SequenceKmerIterable.this.getKmer(pos));
                pos++;
                retVal = true;
            }
            return retVal;
        }

        @Override
        public void forEachRemaining(Consumer<? super String> action) {
            while (pos < limit) {
                action.accept(SequenceKmerIterable.this.getKmer(pos));
                pos++;
            }
        }

        @Override
        public Spliterator<String> trySplit() {
            Spliterator<String> retVal;
            // Note we use a bit-shift for performance.
            int middle = (pos + limit) >> 1;
            if (middle <= pos)
                retVal = null;
            else {
                // Here we have room to split.  "middle" is greater than "pos" and strictly less than "limit".
                // The new splitter ends before middle, and this one moves to it.
                retVal = SequenceKmerIterable.this.new Splitter(pos, middle);
                this.pos = middle;
                // Insure the size stays accurate.
                this.start = middle;
            }
            return retVal;
        }

        @Override
        public long estimateSize() {
            return this.limit - this.start;
        }

        @Override
        public int characteristics() {
            return Spliterator.IMMUTABLE + Spliterator.NONNULL + Spliterator.SIZED + Spliterator.SUBSIZED;
        }

    }

    /**
     * Construct a kmer iterable from a string.
     *
     * @param seq	sequence string to iterate
     * @param k		kmer size
     */
    public SequenceKmerIterable(String seq, int k) {
        this.seqString = seq;
        this.kmerSize = k;
    }

    /**
     * Construct a kmer iterable from a sequence
     *
     * @param seq	sequence to iterate
     * @param k		kmer size
     */
    public SequenceKmerIterable(Sequence seq, int k) {
        this.seqString = seq.getSequence();
        this.kmerSize = k;
    }

    /**
     * @return the position beyond the last permissible kmer start for this sequence
     */
    public int getLimit() {
        return this.seqString.length() - this.kmerSize + 1;
    }

    @Override
    public Iterator<String> iterator() {
        return this.new Iter();
    }

    @Override
    public Spliterator<String> spliterator() {
        return this.new Splitter(0, this.getLimit());
    }

    /**
     * @return the kmer at the specified position.
     *
     * @param pos	starting position of kmer
     */
    public String getKmer(int pos) {
        return this.seqString.substring(pos, pos + this.kmerSize);
    }

    /**
     * @return a stream for these kmers
     */
    public Stream<String> stream() {
        return StreamSupport.stream(this.spliterator(), false);
    }

    /**
     * @return a stream for these kmers
     */

    /**
     * @return a normal or parallel stream for these kmers
     *
     * @param para	TRUE for a parallel stream, else FALSE
     */
    public Stream<String> stream(boolean para) {
        return StreamSupport.stream(this.spliterator(), para);
    }

    /**
     * @return a parallel stream for these kmers
     */
    public Stream<String> parallelStream() {
        return StreamSupport.stream(this.spliterator(), true);

    }

}
