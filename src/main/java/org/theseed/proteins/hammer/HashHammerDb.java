/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
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
    /** map of hammers to genome IDs */
    private Map<String, String> hammerMap;

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

        @Override
        public void updateHammerMap(String fid, String hammer) {
            HashHammerDb.this.hammerMap.put(hammer, fid);
        }

        @Override
        public void createEmptyMap(File inFile) {
            int estimate = (int) (inFile.length() / 30);
            log.info("Estimated hash size for {} is {}.", inFile, estimate);
            HashHammerDb.this.hammerMap = new HashMap<String, String>(estimate);
        }

        @Override
        public void close() {
            // No resources to release
        }

    }

    @Override
    protected void findClosestInternal(CountMap<String> map, Collection<Sequence> seqs, final int kSize) {
        Iterable<String> kIter = KmerSeries.init(seqs, kSize);
        for (String kmer : kIter) {
            String fid = this.hammerMap.get(kmer);
            if (fid != null)
                map.count(Feature.genomeOf(fid));
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
                String fid = this.hammerMap.get(kmer);
                if (fid != null) {
                    // Here we have a hammer hit.  Form the hit descriptor.
                    var hit = new HammerDb.Hit(contigId, len, i, dir, fid, kSize);
                    collection.add(hit);
                }
            }
        }
    }

}
