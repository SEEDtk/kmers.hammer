/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

import org.theseed.counters.WeightMap;
import org.theseed.sequence.KmerSeries;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ParseFailureException;

/**
 * This object implements a thor-hammer database using an in-memory hash.  It is much more memory-intensive
 * than the SQL version, but considerably faster.
 *
 * Each hammer is encoded as a long integer.  The hash is a two-level array, with the first level being a fixed size.
 * If K is the kmer size, then the first level has 4^(K - 15) entries, with a minimum of 1.  For a standard hammer of
 * 20 characters, this is 1024.  Each lower level is an instance of the SubHash class.  The SubHash starts with 31 entries,
 * and each time it fills, it expands to twice the size plus 1 (so that the size is always odd).
 *
 * It is important to note that hammers are never deleted.  This saves us a lot of extra work.
 *
 * @author Bruce Parrello
 *
 */
public class HashHammerDb extends HammerDb {

    // FIELDS
    /** map of hammers to hammer source data */
    private SubHash[] hammerMap;
    /** name of hammer load file */
    private File dbFile;
    /** map of genome IDs to hammer lists */
    private Map<String, HammerArray> genomeMap;
    /** maximum hash size */
    public static final int MAX_HASH = 0x10000000;
    /** starting subhash size */
    public static final int SUBHASH_START = 31;

    /**
     * Construct a hammer database from an input file.  The file must be tab-delimited, with headers,
     * and should contain the hammer sequences in the first column and the feature IDs in the second.
     *
     * @param inFile	input file containing the hammer database
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    public HashHammerDb(File inFile) throws IOException, ParseFailureException {
        this.load(inFile);
        this.dbFile = inFile;
    }

    /**
     * Construct a hammer database for a command processor.  We basically ask for the input file and then
     * load.
     */
    public HashHammerDb(IParms processor) throws ParseFailureException, IOException {
        this.dbFile = processor.getDbFile();
        if (dbFile == null)
            throw new ParseFailureException("File must be specified for in-memory hammer database.");
        this.load(this.dbFile);
    }

    /**
     * This is the subclass for the hash node.  It includes the hammer itself as well as the hammer source.
     */
    protected class SourceNode extends Source {

        /** encoded hammer */
        private long hammerCode;
        /** next source node with the same hash value */
        private SourceNode nextNode;

        /**
         * @param hCode	the encoded hammer in question
         * @param fid		the source feature ID
         * @param str		the hammer strength (varies from 0 to 1)
         */
        protected SourceNode(long hCode, String fid, double str) {
            super(fid, str);
            this.hammerCode = hCode;
            this.nextNode = null;
        }

        /**
         * @return the hammer for this node
         */
        protected String getHammer() {
            return HashHammerDb.this.decode(this.hammerCode);
        }

        /**
         * @return the hammer code
         */
        public long getHammerCode() {
            return this.hammerCode;
        }

        /**
         * @return the nextNode
         */
        protected SourceNode getNextNode() {
            return this.nextNode;
        }

        /**
         * @param nextNode the nextNode to set
         */
        protected void setNextNode(SourceNode nextNode) {
            this.nextNode = nextNode;
        }


    }

    /**
     * This class represents one of the smaller hash tables.  Each is implemented as an array of source nodes, though
     * each node is itself the head of a list of synonyms.  When the number of source nodes gets larger than 3/4 of
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
        private SourceNode[] map;

        /**
         * Construct an empty sub-hash.
         *
         * @param cap	desired capacity
         */
        protected SubHash(int cap) {
            // Here we count on the fact all the pointers start out NULL.
            this.map = new SourceNode[cap];
            this.size = 0;
            this.setMaxCapacity();
        }

        /**
         * Compute the maximum capacity.  If we can't allow the array to grow any more, we disable the
         * maximum-capacity check entirely.
         */
        private void setMaxCapacity() {
            if (this.map.length >= MAX_HASH)
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
            SourceNode[] oldMap = this.map;
            // Here we count on the fact all the pointers start out NULL.
            this.setMaxCapacity();
            this.map = new SourceNode[this.map.length * 2 + 1];
            // We need to copy the nodes from the old map.  Note that we save the next-node
            // pointer before we add to the new map, since it is overridden
            for (int i = 0; i < oldMap.length; i++) {
                SourceNode curr = oldMap[i];
                while (curr != null) {
                    SourceNode next = curr.getNextNode();
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
        private void internalAdd(SourceNode node) {
            int idx = this.computeIdx(node.getHammerCode());
            node.setNextNode(this.map[idx]);
            this.map[idx] = node;
        }

        /**
         * @return the hash index for the specified hammer code
         *
         * @param hammerCode	encoded hammer to check
         */
        private int computeIdx(long hammerCode) {
            // Take the last 30 bits of the hammer and modulo them with the array size.
            int retVal = (int) (hammerCode & 0x3FFFFFFF) % this.map.length;
            return retVal;
        }

        /**
         * Add a new hammer to this hash.
         *
         * @param hammerCode	encoded hammer to add
         * @param fid			source feature ID
         * @param str			hammer strength
         *
         * @return the source node added
         */
        protected SourceNode add(long hammerCode, String fid, double str) {
            // Insure there is room for the new hammer.
            if (this.isFull())
                this.expand();
            // Create a node and add it.
            SourceNode retVal = new SourceNode(hammerCode, fid, str);
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
        protected SourceNode get(long hammerCode) {
            int idx = this.computeIdx(hammerCode);
            SourceNode retVal = this.map[idx];
            while (retVal != null && retVal.getHammerCode() != hammerCode)
                retVal = retVal.getNextNode();
            return retVal;
        }

    }



    /**
     * This is the loader subclass.
     */
    protected class Loader implements HammerDb.ILoader {

        @Override
        public void createEmptyMap(File inFile) {
            // We need to compute the size of the subhash array and fill it with empty subhashes.
            int kSize = HashHammerDb.this.getKmerSize();
            if (kSize == 0)
                throw new IllegalStateException("Cannot create hash hammer map with unknown kmer size.");
            int level1Size = HashHammerDb.getLevel1Size(kSize);
            // Now "level1Size" represents the size of an array that has one slot for every possible bit combination
            // beyond the lowest fifteen base pairs.  We need to fill in the empty hashes.
            HashHammerDb.this.hammerMap = new SubHash[level1Size];
            for (int i = 0; i < level1Size; i++)
                HashHammerDb.this.hammerMap[i] = new SubHash(SUBHASH_START);
            // Finally, create the genome map.
            HashHammerDb.this.genomeMap = new HashMap<String, HammerArray>();
        }

        @Override
        public void updateHammerMap(String fid, String hammer, double str) {
            // Encode the hammer.
            long hammerCode = HashHammerDb.this.encode(hammer);
            // Find the sub-hash and add the hammer to it.
            SubHash subHash = HashHammerDb.this.getSubHash(hammerCode);
            SourceNode node = subHash.add(hammerCode, fid, str);
            // Add the hammer to the genome's hammer array.
            HammerArray gHammers = HashHammerDb.this.genomeMap.computeIfAbsent(node.getGenomeId(),
                    x -> new HammerArray(HashHammerDb.this.getKmerSize()));
            gHammers.add(hammerCode);
        }

        @Override
        public void close() {
        }

    }

    /**
     * @return the sub-hash for a specified hammer
     *
     * @param hammerCode	encoded hammer
     */
    private SubHash getSubHash(long hammerCode) {
        // Compute the index of the subhash.
        int subHashIdx = HashHammerDb.getSubHashIndex(hammerCode);
        SubHash subHash = this.hammerMap[subHashIdx];
        return subHash;
    }

    /**
     * @return the sub-hash index for an encoded hammer
     *
     * @param hammerCode	encoded hammer
     */
    public static int getSubHashIndex(long hammerCode) {
        return (int) (hammerCode >>> 30);
    }

    /**
     * @return the required number of entries in the level-one hash for this kmer size
     *
     * @param kSize		kmer size to use
     */
    public static int getLevel1Size(int kSize) {
        int retVal = 1;
        if (kSize > 15)
            retVal <<= 2 * (kSize - 15);
        return retVal;
    }

    @Override
    protected void findClosestInternal(WeightMap map, Collection<Sequence> seqs, final int kSize) {
        Iterable<String> kIter = KmerSeries.init(seqs, kSize);
        for (String kmer : kIter) {
            HammerDb.Source source = this.getSource(kmer);
            if (source != null)
                this.countHit(map, source, 1);
        }
    }

    /**
     * @return the source for a specified hammer, or NULL if the hammer is not in this database
     *
     * @param kmer		kmer representing a potential hammer
     */
    @Override
    public Source getSource(String kmer) {
        Source retVal;
        long hammerCode = this.encode(kmer);
        if (hammerCode < 0) {
            // Here we have an invalid character, so the answer is always no.
            retVal = null;
        } else {
            SubHash subHash = this.getSubHash(hammerCode);
            retVal = subHash.get(hammerCode);
        }
        return retVal;
    }

    @Override
    protected ILoader getLoader() {
        return this.new Loader();
    }

    @Override
    protected void findHitsInternal(Collection<Hit> collection, Collection<Sequence> seqs, int kSize, boolean dir) {
        log.debug("Scanning {} sequences for hammer hits with strand flag {}.", seqs.size(), dir);
        for (Sequence seq : seqs) {
            String dna = seq.getSequence();
            final int len = dna.length();
            String contigId = seq.getLabel();
            final int n = len - kSize;
            // We loop through the sequence with a character index, since we need the location of the hit.
            for (int i = 0; i <= n; i++) {
                String kmer = dna.substring(i, i + kSize);
                HammerDb.Source source = this.getSource(kmer);
                if (source != null) {
                    // Here we have a hammer hit.  Form the hit descriptor.
                    var hit = new HammerDb.Hit(contigId, len, i, dir, source.getFid(), kSize, source.getStrength());
                    collection.add(hit);
                }
            }
        }
    }

    @Override
    public File getLoadFile() {
        return this.dbFile;
    }

    @Override
    protected void findHammersInternal(HashSet<String> hammerSet, String seq, int kSize) {
        SequenceKmerIterable iter = new SequenceKmerIterable(seq, kSize);
        for (String kmer : iter) {
            if (this.getSource(kmer) != null)
                hammerSet.add(kmer);
        }
    }

    @Override
    public Map<String, Source> findGenomeHammers(String genomeId) {
        Map<String, Source> retVal = new HashMap<String, Source>();
        // We keep a list of all the hammers for each genome.  We just need to find it.
        HammerArray hammers = this.genomeMap.get(genomeId);
        if (hammers != null) {
            // We found it, so run through all the hammers.
            for (String hammer : hammers) {
                Source node = this.getSource(hammer);
                retVal.put(hammer, node);
            }
        }
        return retVal;
    }

}
