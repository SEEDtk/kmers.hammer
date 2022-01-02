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
        public void updateHammerMap(String genome, String hammer) {
            HashHammerDb.this.hammerMap.put(hammer, genome);
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
    protected CountMap<String> findClosestInternal(Collection<Sequence> seqs, final int kSize) {
        CountMap<String> retVal = new CountMap<String>();
        for (Sequence seq : seqs) {
            String dna = seq.getSequence();
            final int n = dna.length() - kSize;
            for (int i = 0; i < n; i++) {
                String kmer = dna.substring(i, i + kSize);
                String genomeId = this.hammerMap.get(kmer);
                if (genomeId != null)
                    retVal.count(genomeId);
            }
        }
        return retVal;
    }

    @Override
    protected ILoader getLoader() {
        return this.new Loader();
    }

}
