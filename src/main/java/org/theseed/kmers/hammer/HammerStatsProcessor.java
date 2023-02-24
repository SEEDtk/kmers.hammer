/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseHammerUsageProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This method reads a FASTA file and counts the number of times each hammer was found.  It will output the source feature
 * ID and the number of occurrences for each such hammer found.
 *
 * The positional parameter should be the name of the FASTA file to use as input.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 *
 * --hType		type of hammer database (default MEMORY)
 * --method		voting method to use (default COUNT)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type
 *
 * @author Bruce Parrello
 *
 */
public class HammerStatsProcessor extends BaseHammerUsageProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerStatsProcessor.class);
    /** map of hammers to info counts */
    private Map<String, HammerCount> countMap;

    /**
     * This object is essentially a mutable integer ready for counting.
     */
    protected static class HammerCount {

        /** count value */
        private int count;

        /**
         * Create a blank, empty hammer count.
         */
        protected HammerCount() {
            this.count = 0;
        }

        /**
         * Increment a hammer count.
         */
        protected void count() {
            this.count++;
        }

        /**
         * @return the total count
         */
        protected int getCount() {
            return this.count;
        }

    }

    // COMMAND-LINE OPTIONS

    /** name of the input FASTA file */
    @Argument(index = 0, metaVar = "input.fa", usage = "name of the input FASTA file")
    private File inFile;

    @Override
    protected void setHammerDefaults() {
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input FASTA file " + this.inFile + " is not found or unreadable.");
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        // Create the count map.
        this.countMap = new HashMap<String, HammerCount>(10000);
        // Initialize the input counter.
        int inCount = 0;
        log.info("Processing input.");
        // Open the FASTA input file.
        try (FastaInputStream inStream = new FastaInputStream(this.inFile)) {
            // Loop through the file.
            for (Sequence seq : inStream) {
                inCount++;
                Set<String> hammers = hammerDb.findHammers(seq.getSequence());
                for (String hammer : hammers) {
                    HammerCount hammerCount = this.countMap.computeIfAbsent(hammer, x -> new HammerCount());
                    hammerCount.count();
                }
                if (log.isInfoEnabled() && inCount % 50000 == 0)
                    log.info("{} sequences processed. {} hammers found.", inCount, this.countMap.size());
            }
        }
        // Now we produce the report.
        log.info("Writing report.");
        writer.println("hammer\tfid\tgenome_id\tcount");
        for (var countInfo : this.countMap.entrySet()) {
            String hammer = countInfo.getKey();
            HammerDb.Source source = hammerDb.getSource(hammer);
            writer.println(hammer + "\t" + source.getFid() + "\t" + source.getGenomeId() + "\t" + Integer.toString(countInfo.getValue().getCount()));
        }

    }

}
