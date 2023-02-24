/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Iterator;

/**
 * This implements a variable-length list of hammers.  It looks like a string list, but in fact it is an array of
 * long integers.
 *
 * @author Bruce Parrello
 *
 */
public class HammerArray implements Iterable<String> {

    // FIELD
    /** array of encoded hammers */
    private long[] hammers;
    /** number of hammers in the array */
    private int size;
    /** kmer size for the hammers */
    private int kSize;
    /** default array size */
    private static final int DEFAULT_SIZE = 1000;
    /** array increment size */
    private static final int DEFAULT_INCR = 500;

    /**
     * Create a new, empty hammer array.
     */
    public HammerArray(int kSize) {
        this.hammers = new long[DEFAULT_SIZE];
        this.kSize = kSize;
        this.size = 0;
    }

    /**
     * This class is an iterator for the array.
     */
    protected class Iter implements Iterator<String> {

        /** position in the array or the next value */
        private int idx;

        /**
         * Initialize the iterator.
         */
        private Iter() {
            this.idx = 0;
        }

        @Override
        public boolean hasNext() {
            return this.idx < HammerArray.this.size;
        }

        @Override
        public String next() {
            String retVal = HammerMap.decode(HammerArray.this.hammers[this.idx], HammerArray.this.kSize);
            this.idx++;
            return retVal;
        }

    }

    /**
     * @return the size of the array
     */
    public int size() {
        return this.size;
    }

    /**
     * @return an iterator through the hammers
     */
    @Override
    public Iterator<String> iterator() {
        return this.new Iter();
    }

    /**
     * Add a new hammer to the array.
     *
     * @param hammerCode	encoded hammer to add
     */
    protected void add(long hammerCode) {
        if (this.size >= this.hammers.length) {
            // Here we need to grow the array.
            long[] newArray = new long[this.size + DEFAULT_INCR];
            System.arraycopy(this.hammers, 0, newArray, 0, this.size);
            this.hammers = newArray;
        }
        this.hammers[this.size] = hammerCode;
        this.size++;
    }

   /**
    * @return the hammer at the specified position in the array
    *
    * @param index	position to check
    */
   public String get(int index) {
       if (index < 0 || index >= this.size)
           throw new IndexOutOfBoundsException("Attempt to access position " + Integer.toString(index) + " in hammer array with "
                   + Integer.toString(this.size) + " values.");
       String retVal = HammerMap.decode(this.hammers[index], this.kSize);
       return retVal;
    }

}
