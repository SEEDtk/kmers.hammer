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

import org.theseed.basic.ParseFailureException;
import org.theseed.sequence.ISequence;
import org.theseed.sequence.KmerSeries;
import org.theseed.sequence.Sequence;

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
    private HammerMap<Source> hammerMap;
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
     * This is the loader subclass.
     */
    protected class Loader implements HammerDb.ILoader {

        @Override
        public void createEmptyMap(File inFile) {
            // We need to compute the size of the subhash array and fill it with empty subhashes.
            int kSize = HashHammerDb.this.getKmerSize();
            if (kSize == 0)
                throw new IllegalStateException("Cannot create hash hammer map with unknown kmer size.");
            // Create the main hammer map.
            HashHammerDb.this.hammerMap = new HammerMap<Source>(kSize);
            // Finally, create the genome map.
            HashHammerDb.this.genomeMap = new HashMap<String, HammerArray>();
        }

        @Override
        public void updateHammerMap(String fid, String roleId, String hammer, double str) {
            // Add the hammer to the main map.
            Source node = new Source(fid, roleId, str);
            HashHammerDb.this.hammerMap.put(hammer, node);
            // Add the hammer to the genome's hammer array.
            HammerArray gHammers = HashHammerDb.this.genomeMap.computeIfAbsent(node.getGenomeId(),
                    x -> new HammerArray(HashHammerDb.this.getKmerSize()));
            gHammers.add(HashHammerDb.this.hammerMap.encode(hammer));
        }

        @Override
        public void close() {
            log.info("Hash load factor is {}. Overload factor is {}.", HashHammerDb.this.hammerMap.loadFactor(),
                    HashHammerDb.this.hammerMap.overloadFactor());
        }

    }

   @Override
    protected void findClosestInternal(ScoreMap map, Collection<Sequence> seqs, final int kSize) {
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
        return this.hammerMap.get(kmer);
    }

    @Override
    protected ILoader getLoader() {
        return this.new Loader();
    }

    @Override
    protected void findHitsInternal(Collection<Hit> collection, Collection<? extends ISequence> seqs, int kSize, boolean dir) {
        log.debug("Scanning {} sequences for hammer hits with strand flag {}.", seqs.size(), dir);
        for (ISequence seq : seqs) {
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
                    var hit = new HammerDb.Hit(contigId, len, i, dir, source.getFid(), source.getRole(), kSize, source.getStrength());
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
