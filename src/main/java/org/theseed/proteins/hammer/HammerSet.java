/**
 *
 */
package org.theseed.proteins.hammer;

/**
 * This object manages a set of hammers.  The set looks like a list of strings, but in fact everything is
 * encoded as a long integer.  It uses the same basic two-level structure as HashHammerDb, but the content
 * is different.
 *
 * @author Bruce Parrello
 *
 */
public class HammerSet {

    // FIELDS
    /** array of sub-hashes */
    private SubHash[] hammerMap;
    /** hammer kmer size */
    private int kmerSize;

    /**
     * This object represents a single node in one of the smaller hash tables.
     */
    protected static class HammerNode {

        /** code for this hammer */
        private long hammerCode;
        /** next node with the same hash value */
        private HammerNode nextNode;

        /**
         * Construct a node for an encoded hammer
         *
         * @param code		hammer code
         */
        public HammerNode(long code) {
            this.hammerCode = code;
        }
        /**
         * @return the next synonom for this hammer, or NULL if none
         */
        public HammerNode getNextNode() {
            return this.nextNode;
        }
        /**
         * @return the hammer code
         */
        public long getHammerCode() {
            return this.hammerCode;
        }
        /**
         * Specify the next synonym node for this hammer hash code.
         *
         * @param hammerNode	next synonym node
         */
        public void setNextNode(HammerNode hammerNode) {
            this.nextNode = hammerNode;
        }

    }

    /**
     * This class represents one of the smaller hash tables.  Each is implemented as an array of hammer nodes, though
     * each node is itself the head of a list of synonyms.  When the number of hammer nodes gets larger than 3/4 of
     * the array size, then if the array size is less than the maximum, the array is made bigger.
     * The hash function is based on the last 30 bits of the hammer code.
     *
     * @author Bruce Parrello
     *
     */
    protected class SubHash {

        /** maximum acceptable capacity */
        private int maxCapacity;
        /** current number of hammers in this table */
        private int size;
        /** array of source nodes */
        private HammerNode[] map;

        /**
         * Construct an empty sub-hash.
         *
         * @param cap	desired capacity
         */
        protected SubHash(int cap) {
            // Here we count on the fact all the pointers start out NULL.
            this.map = new HammerNode[cap];
            this.size = 0;
            this.setMaxCapacity();
        }

        /**
         * Compute the maximum capacity.  If we can't allow the array to grow any more, we disable the
         * maximum-capacity check entirely.
         */
        private void setMaxCapacity() {
            if (this.map.length >= HashHammerDb.MAX_HASH)
                this.maxCapacity = 0;
            else
                this.maxCapacity = (this.map.length * 3) / 4;
        }

        /**
         * @return TRUE if there is no room for another node in this sub-hash
         */
        protected boolean isFull() {
            return (this.maxCapacity > 0 && this.size >= this.maxCapacity);
        }

        /**
         * Expand this sub-hash to make more room.
         */
        private void expand() {
            HammerNode[] oldMap = this.map;
            // Here we count on the fact all the pointers start out NULL.
            this.setMaxCapacity();
            this.map = new HammerNode[this.map.length * 2 + 1];
            // We need to copy the nodes from the old map.  Note that we save the next-node
            // pointer before we add to the new map, since it is overridden
            for (int i = 0; i < oldMap.length; i++) {
                HammerNode curr = oldMap[i];
                while (curr != null) {
                    HammerNode next = curr.getNextNode();
                    this.internalAdd(curr);
                    curr = next;
                }
            }
        }

        /**
         * Add a new node to this hash.  This is an internal-use-only method.   We don't update the size,
         * since it is used during splits as well as during random adds.
         *
         * @param node		node to add
         */
        private void internalAdd(HammerNode node) {
            int idx = this.computeIdx(node.getHammerCode());
            node.setNextNode(this.map[idx]);
            this.map[idx] = node;
        }

        /**
         * @return the hash index for the specified hammer code
         *
         * @param hammerCode	encoded hammer to check
         */
        protected int computeIdx(long hammerCode) {
            // Take the last 30 bits of the hammer and modulo them with the array size.
            int retVal = (int) (hammerCode & 0x3FFFFFFF) % this.map.length;
            return retVal;
        }

        /**
         * Add a new hammer to this hash.
         *
         * @param hammerCode	encoded hammer to add
         *
         * @return the source node added
         */
        protected HammerNode add(long hammerCode) {
            // Insure there is room for the new hammer.
            if (this.isFull())
                this.expand();
            // Create a node and add it.
            HammerNode retVal = new HammerNode(hammerCode);
            this.internalAdd(retVal);
            this.size++;
            return retVal;
        }

        /**
         * Find a hammer in this hash.
         *
         * @param hammerCode	encoded hammer to find
         *
         * @return the source node for the hammer, or NULL if there was none
         */
        protected HammerNode get(long hammerCode) {
            int idx = this.computeIdx(hammerCode);
            HammerNode retVal = this.map[idx];
            while (retVal != null && retVal.getHammerCode() != hammerCode)
                retVal = retVal.getNextNode();
            return retVal;
        }

    }

    /**
     * Construct an empty hammer set.
     *
     * @param	kSize	hammer kmer size
     */
    public HammerSet(int kSize) {
        this.kmerSize = kSize;
        int level1Size = HashHammerDb.getLevel1Size(kSize);
        this.hammerMap = new SubHash[level1Size];
        for (int i = 0; i < level1Size; i++)
            this.hammerMap[i] = new SubHash(HashHammerDb.SUBHASH_START);
    }

    /**
     * Add a hammer to the set.
     *
     * @param hammer	hammer to add (un-encoded)
     *
     * @return TRUE if we added the hammer, FALSE if it was already present
     */
    public boolean add(String hammer) {
        long code = HammerMap.encode(hammer, this.kmerSize);
        SubHash subHash = this.hammerMap[HashHammerDb.getSubHashIndex(code)];
        // Find the hammer in this sub-hash
        HammerNode node = subHash.get(code);
        // If we didn't find it, add it.
        boolean retVal = (node == null);
        if (retVal)
            subHash.add(code);
        return retVal;
    }

    /**
     * Check to determine if a hammer is in this set.
     *
     * @param hammer	hammer to check (un-encoded)
     *
     * @return TRUE if the hammer is in the set, else FALSE
     */
    public boolean contains(String hammer) {
    	boolean retVal;
        long code = HammerMap.encode(hammer, this.kmerSize);
        if (code < 0)
        	retVal = false;
        else {
	        SubHash subHash = this.hammerMap[HashHammerDb.getSubHashIndex(code)];
	        // Find the hammer in this sub-hash
	        HammerNode node = subHash.get(code);
	        retVal = (node != null);
        }
        return retVal;
    }


}
