/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.NoSuchElementException;
import java.util.function.Consumer;
import java.util.function.Function;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This object implements a hammer map.  Because it can hold more than two billion entries, it is
 * not a java map.  Adds and removes are thread-safe, but traversals are not.
 *
 * The hash is implemented in two levels, and all the hammers are stored encoded.  The main hash is a
 * simple lookup array based on the first few bits of the encoded hash.  Each entry in the array is a
 * SubHash with chaining.
 *
 * @author Bruce Parrello
 *
 */
public class HammerMap<T> implements Iterable<Map.Entry<String, T>> {


    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerMap.class);
   /** hammer kmer size */
    private int kmerSize;
    /** fixed array of sub-hashes */
    private ArrayList<SubHash> hammerMap;
    /** conversion array for encoding/decoding */
    protected static final char[] CONVERTER = new char[] { 'a', 'c', 'g', 't' };
    /** number of bits occupied by the second-level part of the hammer code */
    protected static final int LOWER_BITS = (Integer.SIZE - 1) / 2 * 2;
    /** the number of base pairs occupied by the second-level part of the hammer code */
    protected static final int LOWER_CHARS = LOWER_BITS / 2;
    /** bit mask for the second-level part of the hammer code */
    protected static final int LOWER_BIT_MASK = (1 << LOWER_BITS) - 1;
    /** load factor for sub-hashes */
    protected static final double LOAD_FACTOR = 0.75;
    /** initial size for sub-hashes */
    protected static final int START_SIZE = 31;
    /** initial maximum capacity for sub-hashes */
    protected static final int START_MAX = (int) (START_SIZE * LOAD_FACTOR);


    protected class Node implements Map.Entry<String, T> {

        /** hammer code for this node */
        private long hammerCode;
        /** value of this node */
        private T value;
        /** next node in this synonym chain */
        private Node nextNode;

        protected Node(long code, T item) {
            this.hammerCode = code;
            this.value = item;
            this.nextNode = null;
        }

        /**
         * @return TRUE if the specified hammer code is in this node, else FALSE
         *
         * @param code	hammer code to check
         */
        protected boolean matches(long code) {
            return (code == this.hammerCode);
        }

        @Override
        public String getKey() {
            return HammerMap.this.decode(this.hammerCode);
        }

        @Override
        public T getValue() {
            return this.value;
        }

        @Override
        public T setValue(T newValue) {
            this.value = newValue;
            return this.value;
        }

        /**
         * Specify the next node in this one's synonym chain.
         *
         * @param node		node to place after this one
         */
        private void setNextNode(Node node) {
            this.nextNode = node;
        }

        /**
         * @return the next node in this one's synonym chain, or NULL if there is none
         */
        private Node getNextNode() {
            return this.nextNode;
        }

    }


    protected class SubHash {

        /** number of hammers in this sub-hash */
        private int size;
        /** desired maximum capacity for this sub-hash, or 0 if it should be allowed to grow unlimited */
        private int maxCapacity;
        /** array of synonym chains; we use an object array because Node is generic */
        private Object[] table;

        /**
         * Construct a new, empty subhash.
         */
        protected SubHash() {
            this.size = 0;
            this.maxCapacity = START_MAX;
            this.table = new Object[START_SIZE];
        }

        /**
         * @return the number of hammers in this sub-hash
         */
        protected int size() {
            return this.size;
        }

        /**
         * @return the number of synonym chains in this sub-hash
         */
        private int capacity() {
            return this.table.length;
        }

        /**
         * Search a synonym chain for a matching node.  Return NULL if there is none.
         *
         * @param code			encoded hammer to find
         * @param subHash		sub-hash containing the desired
         * @param subIdx		index of the hammer's synonym chain
         *
         * @return the node found, or NULL if it is not present
         */
        private Node findMatch(long code, int subIdx) {
            Node retVal = this.getChain(subIdx);
            while (retVal != null && ! retVal.matches(code))
                retVal = retVal.getNextNode();
            return retVal;
        }

        /**
         * @return the head of the synonym chain at the specified index, or NULL if the chain is empty
         *
         * @param subIdx	index of desired chain
         */
        @SuppressWarnings("unchecked")
        private Node getChain(int subIdx) {
            return (Node) this.table[subIdx];
        }

        /**
         * InterNal-use only method for adding a hammer/value pair when the chain and the
         * array index are known.
         *
         * @param code		hammer code
         * @param subIdx	index of the synonym chain for this hammer
         * @param value		value to assign to the hammer
         */
        private void addInternal(long code, int subIdx, T value) {
            // Insure there is room in this sub-hash.
            if (this.maxCapacity > 0 && this.size >= this.maxCapacity) {
                this.expand();
                subIdx = this.getSubIdx(code);
            }
            // Create the node and add it.
            Node newNode = new Node(code, value);
            this.addNode(subIdx, newNode);
        }

        /**
         * This adds the hammer-value node to a synonym chain.
         *
         * @param subIdx	index of the synonym chain for the node
         * @param newNode	node to add
         */
        private void addNode(int subIdx, Node newNode) {
            Node newNext = this.getChain(subIdx);
            newNode.setNextNode(newNext);
            this.table[subIdx] = newNode;
            this.size++;
        }

        /**
         * Make this sub-hash bigger so that a new node can fit.
         */
        private void expand() {
            // Compute the new maximum capacity.
            long newSize = ((long) this.capacity()) * 2 + 1;
            if (newSize > Integer.MAX_VALUE) {
                // We can't expand.  Insure we stop trying.
                this.maxCapacity = 0;
            } else {
                // Create the bigger array.
                var oldTable = this.table;
                this.table = new Object[(int) newSize];
                this.maxCapacity = (int) (newSize * LOAD_FACTOR);
                // Reset ourselves to empty, as we are going to be re-inserting everything.
                this.size = 0;
                // Loop through the table, copying synonym chains.
                for (int i = 0; i < oldTable.length; i++) {
                    @SuppressWarnings("unchecked")
                    Node oldNode = (Node) oldTable[i];
                    // Loop through the synonym chain, copying to the new table.
                    while (oldNode != null) {
                        int subIdx = this.getSubIdx(oldNode.hammerCode);
                        Node nextNode = oldNode.getNextNode();
                        this.addNode(subIdx, oldNode);
                        oldNode = nextNode;
                    }
                }
            }
        }

        /**
         * @return the table entry relevant to the specified hammer code
         *
         * @param code	encoded hammer
         */
        private int getSubIdx(long code) {
            return (int) (code & LOWER_BIT_MASK) % this.capacity();
        }

        /**
         * Remove the specified hammer from this hash, and return its value.
         *
         * @param code	encoded hammer to remove
         *
         * @return the value of the hammer, or NULL if it is not found
         */
        private T remove(long code) {
            T retVal = null;
            int subIdx = this.getSubIdx(code);
            Node curr = this.getChain(subIdx);
            if (curr.matches(code)) {
                // The first node is us (the most common case), so we simply replace it
                // with its successor.
                this.table[subIdx] = curr.getNextNode();
                retVal = curr.getValue();
            } else {
                // Search for this hammer in the chain.
                Node prev = curr;
                curr = curr.nextNode;
                while (curr != null && ! curr.matches(code)) {
                    prev = curr;
                    curr = curr.getNextNode();
                }
                if (curr != null) {
                    // Here we found a match.
                    prev.setNextNode(curr.getNextNode());
                    retVal = curr.getValue();
                }
            }
            return retVal;
        }

    }

    /**
     * This class implements an iterator through the hash.  Terrible things will happen
     * if the hash is modified while iterating.
     */
    public class Iter implements Iterator<Entry<String, T>> {

        /** current sub-hash */
        private SubHash subHash;
        /** index of next sub-hash */
        private int mainIdx;
        /** index of next table entry in sub-hash */
        private int subIdx;
        /** next node to return */
        private Node next;

        /**
         * Initialize this iterator.
         */
        public Iter() {
            // Position on the first chain in the first sub-hash.
            this.subHash = HammerMap.this.hammerMap.get(0);
            this.subIdx = 0;
            // Denote the next sub-hash will be the second one.
            this.mainIdx = 1;
            // Denote we do not have a next-node yet.
            this.next = null;
            // Find the first nonempty chain.
            this.findNextChain();
        }

        /**
         * Position on the first node in the next non-empty chain
         * in the map.
         */
        private void findNextChain() {
            final int n = HammerMap.this.hammerMap.size();
            // Scan the current subhash for a nonempty chain.
            this.scanSubHash();
            // If we did not find anything, loop until we do find
            // something or we run out of subhashes.
            while (this.next == null && this.mainIdx < n) {
                // Get the next sub-hash.
                this.subHash = HammerMap.this.hammerMap.get(this.mainIdx);
                this.mainIdx++;
                this.subIdx = 0;
                // Scan it for a nonempty chain.
                this.scanSubHash();
            }
        }

        /**
         * Scan the current subhash for the next non-empty chain.
         */
        private void scanSubHash() {
            final int n = this.subHash.capacity();
            while (this.subIdx < n && this.next == null) {
                this.next = this.subHash.getChain(this.subIdx);
                this.subIdx++;
            }
        }

        @Override
        public boolean hasNext() {
            return this.next != null;
        }

        @Override
        public Entry<String, T> next() {
            if (this.next == null)
                throw new NoSuchElementException("Attempt to read past end of hammer map.");
            var retVal = this.next;
            // Read ahead to the next node.  This tells us if there is one, and preps us
            // for our next iteration.
            this.next = retVal.getNextNode();
            if (this.next == null)
                this.findNextChain();
            return retVal;
        }

    }


    /**
     * Encode a hammer into a long integer.  The hammer string should match the kmer size.  If it contains
     * any ambiguity characters, it is treated as invalid.
     *
     * @param hammer	hammer to encode
     * @param kSize		kmer size for the hammer
     *
     * @return a long integer representation of the hammer, or a negative value if the hammer is invalid
     */
    public static long encode(String hammer, final int kSize) {
        long retVal = 0;
        if (hammer.length() != kSize)
            retVal = -1;
        else for (int i = 0; i < kSize && retVal >= 0; i++) {
            retVal <<= 2;
            switch (hammer.charAt(i)) {
            case 'A' :
            case 'a' :
                break;
            case 'C' :
            case 'c' :
                retVal |= 1;
                break;
            case 'G' :
            case 'g' :
                retVal |= 2;
                break;
            case 'T' :
            case 't' :
                retVal |= 3;
                break;
            default:
                retVal = -1;
            }
        }
        return retVal;
    }

    /**
     * Decode a hammer from a long integer.
     *
     * @param coded		long integer to decode
     * @param kSize		kmer size of the hammer
     *
     * @return a string representation of the hammer, or NULL if the hammer code is invalid
     */
    public static String decode(long coded, final int kSize) {
        String retVal;
        if (coded < 0)
            retVal = null;
        else {
            char[] buffer = new char[kSize];
            for (int i = kSize - 1; i >= 0; i--) {
                buffer[i] = CONVERTER[(int) (coded & 3)];
                coded >>= 2;
            }
            retVal = String.valueOf(buffer);
        }
        return retVal;
    }

    /**
     * @return the number of occurrences of the most common base pair in a sequence
     */
    public static int commonBaseCount(String seq) {
        // Here we count on the fact the array is pre-filled with zeros.
        int[] counts = new int[4];
        for (int i = 0; i < seq.length(); i++) {
            switch (seq.charAt(i)) {
            case 'A' :
            case 'a' :
                counts[0]++;
                break;
            case 'C' :
            case 'c' :
                counts[1]++;
                break;
            case 'G' :
            case 'g' :
                counts[2]++;
                break;
            case 'T' :
            case 't' :
                counts[3]++;
                break;
            }
        }
        int retVal = counts[0];
        for (int i = 1; i < 4; i++)
            if (counts[i] > retVal) retVal = counts[i];
        return retVal;
    }

    /**
     * Construct an empty hammer map.
     *
     * @param kmerSize	the length of a hammer, in base pairs
     */
    public HammerMap(int kmerSize) {
        // Save the kmer size.
        this.kmerSize = kmerSize;
        // Determine how big a level-1 array we need.
        int highLetters = kmerSize - LOWER_CHARS;
        int level1Size = (highLetters <= 0 ? 1 : 1 << (highLetters * 2));
        // Create the level-1 array.
        this.hammerMap = new ArrayList<SubHash>(level1Size);
        for (int i = 0; i < level1Size; i++)
            this.hammerMap.add(new SubHash());
    }

    /**
     * @return the code for a hammer from this map
     *
     * @param hammer	hammer to encode
     */
    public long encode(String hammer) {
        return encode(hammer, this.kmerSize);
    }

    /**
     * @return the string representation of a hammer code from this map, or NULL if the hammer code is invalid
     *
     * @param hammerCode	hammer code
     */
    public String decode(long hammerCode) {
        return decode(hammerCode, this.kmerSize);
    }

    /**
     * @return the object associated with a specified hammer, or NULL if it is not present.
     *
     * @param hammer		hammer to check
     */
    public T get(String hammer) {
        T retVal;
        Node node = this.get(this.encode(hammer));
        if (node == null)
            retVal = null;
        else
            retVal = node.getValue();
        return retVal;
    }

    /**
     * @return the node for a specified hammer, or NULL if it is not in the map
     *
     * @param hammerCode	encoded hammer
     */
    private Node get(long encode) {
        Node retVal;
        if (encode < 0)
            retVal = null;
        else {
            SubHash subHash = this.getSubHash(encode);
            int subIdx = subHash.getSubIdx(encode);
            retVal = subHash.findMatch(encode, subIdx);
        }
        return retVal;
    }

    /**
     * Add a hammer to this map.
     *
     * @param hammer	hammer to add
     * @param value		value to associate with the hammer
     *
     * @return TRUE if the hammer was added, FALSE if the hammer's value was replaced
     */
    public boolean put(String hammer, T value) {
        boolean retVal;
        long code = this.encode(hammer);
        if (code < 0)
            throw new IllegalArgumentException("Invalid hammer \"" + hammer + "\" specified.");
        SubHash subHash = this.getSubHash(code);
        // Here is where we make it all thread-safe.  We cannot do ANYTHING to the sub-hash unless
        // it is locked to us.
        synchronized (subHash) {
            int subIdx = subHash.getSubIdx(code);
            Node nodeFound = subHash.findMatch(code, subIdx);
            if (nodeFound == null) {
                subHash.addInternal(code, subIdx, value);
                retVal = true;
            } else {
                nodeFound.setValue(value);
                retVal = false;
            }
        }
        return retVal;
    }

    /**
     * Insert a hammer and its value in the map, or update the value if it is already present.
     *
     * @param hammer		hammer of interest
     * @param oldFunction	operation to update the value if it is found (NULL to do nothing)
     * @param newFunction	function to create the value if it is not found
     *
     * @return TRUE if we created the value, FALSE if we simply updated it
     */
    public boolean update(String hammer, Consumer<T> oldConsumer, Function<String, ? extends T> newFunction) {
        boolean retVal;
        long code = this.encode(hammer);
        if (code < 0)
            throw new IllegalArgumentException("Invalid hammer \"" + hammer + "\" specified.");
        else {
            SubHash subHash = this.getSubHash(code);
            // Here is where we make it all thread-safe.  We cannot do ANYTHING to the sub-hash unless
            // it is locked to us.
            synchronized (subHash) {
                int subIdx = subHash.getSubIdx(code);
                Node nodeFound = subHash.findMatch(code, subIdx);
                if (nodeFound == null) {
                    T newValue = newFunction.apply(hammer);
                    subHash.addInternal(code, subIdx, newValue);
                    retVal = true;
                } else {
                    retVal = false;
                    if (oldConsumer != null)
                        oldConsumer.accept(nodeFound.getValue());
                }
            }
        }
        return retVal;
    }

    /**
     * @return the sub-hash for the specified hammer code.
     *
     * @param code	encoded hammer (MUST be >= 0)
     */
    private SubHash getSubHash(long code) {
        int idx = (int) (code >>> LOWER_BITS);
        return this.hammerMap.get(idx);
    }

    /**
     * Remove a hammer from the hash.
     *
     * @param hammer	hammer to be removed
     *
     * @return the value of the hammer, or NULL if the hammer was not in the hash
     */
    public T remove(String hammer) {
        T retVal = null;
        long code = this.encode(hammer);
        if (code >= 0) {
            SubHash subHash = this.getSubHash(code);
            synchronized (subHash) {
                retVal = subHash.remove(code);
            }
        }
        return retVal;
    }

    /**
     * @return the number of hammers in the map
     */
    public long size() {
        long retVal = 0;
        for (SubHash subHash : this.hammerMap)
            retVal += subHash.size();
        return retVal;
    }

    /**
     * @return the mean length of an occupied chain
     */
    public double overloadFactor() {
        long totalLength = 0;
        int chains = 0;
        for (SubHash subHash : this.hammerMap) {
            final int n = subHash.capacity();
            for (int i = 0; i < n; i++) {
                Node node = subHash.getChain(i);
                if (node != null) {
                    chains++;
                    while (node != null) {
                        totalLength++;
                        node = node.getNextNode();
                    }
                }
            }
        }
        double retVal = 0.0;
        if (chains > 0)
            retVal = ((double) totalLength) / chains;
        return retVal;
    }

    /**
     * @return the ratio of occupied chains to chain slots
     */
    public double loadFactor() {
        int chains = 0;
        long buckets = 0;
        for (SubHash subHash : this.hammerMap) {
            final int n = subHash.capacity();
            for (int i = 0; i < n; i++) {
                Node node = subHash.getChain(i);
                buckets++;
                if (node != null)
                    chains++;
            }
        }
        double retVal = 0.0;
        if (buckets > 0)
            retVal = chains / (double) buckets;
        return retVal;
    }


    @Override
    public Iterator<Map.Entry<String, T>> iterator() {
        return this.new Iter();
    }

}
