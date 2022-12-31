/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.theseed.counters.WeightMap;
import org.theseed.sequence.KmerSeries;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ParseFailureException;

/**
 * This object implements a thor-hammer database using an in-memory hash.  It should only be used
 * for small ones.
 *
 * @author Bruce Parrello
 *
 */
public class HashHammerDb extends HammerDb {

    // FIELDS
    /** map of hammers to hammer source data */
    private Map<String, Source> hammerMap;
    /** maximum hash size */
    private static final int MAX_HASH = 0x10000000;

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
    }

    public HashHammerDb(IParms processor) throws ParseFailureException, IOException {
        File dbFile = processor.getDbFile();
        if (dbFile == null)
            throw new ParseFailureException("File must be specified for in-memory hammer database.");
        this.load(processor.getDbFile());
    }

    /**
     * This is the loader subclass.
     */
    protected class Loader implements HammerDb.ILoader {

        /** estimated hash size */
        private int estimate;

        @Override
        public void updateHammerMap(String fid, String hammer, double str) {
            HashHammerDb.this.hammerMap.put(hammer, new HammerDb.Source(fid, str));
        }

        @Override
        public void createEmptyMap(File inFile) {
            // We do some fancy math here to insure we have a working hash even if the
            // number of hammers is insane.  We cap the estimated size around 2^28, and
            // if the hash is going to be big, we set the load factor to 1 to prevent
            // splitting.
            float loadFactor = 0.75f;
            this.estimate = (int) (inFile.length() / 50);
            if (this.estimate <= 0 || this.estimate >= MAX_HASH) {
                this.estimate = MAX_HASH;
                loadFactor = 1.00f;
            }
            log.info("Estimated hash size for {} is {}.", inFile, this.estimate);
            HashHammerDb.this.hammerMap = new HashMap<String, HammerDb.Source>(estimate, loadFactor);
        }

        @Override
        public void close() {
            if (log.isInfoEnabled()) {
                double hammerCount = HashHammerDb.this.hammerMap.size();
                log.info("Final hash metrics are {} hammers with {} capacity.", hammerCount,
                    this.estimate / hammerCount);
            }
        }

    }

    @Override
    protected void findClosestInternal(WeightMap map, Collection<Sequence> seqs, final int kSize) {
        Iterable<String> kIter = KmerSeries.init(seqs, kSize);
        for (String kmer : kIter) {
            HammerDb.Source source = this.hammerMap.get(kmer);
            if (source != null)
                this.countHit(map, source);
        }
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
                HammerDb.Source source = this.hammerMap.get(kmer);
                if (source != null) {
                    // Here we have a hammer hit.  Form the hit descriptor.
                    var hit = new HammerDb.Hit(contigId, len, i, dir, source.getFid(), kSize, source.getStrength());
                    collection.add(hit);
                }
            }
        }
    }

}
