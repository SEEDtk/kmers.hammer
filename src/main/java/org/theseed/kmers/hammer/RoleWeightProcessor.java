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

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.io.TabbedLineReader.Line;

/**
 * This command reads a hammer file and changes the weights to scale them by the number of hammers per role/genome combination.  So, if
 * a genome has 100 hammers for a role, the weight would be BASE / 100, where BASE is the average number of hammers per role/genome
 * pair.  This is a two-pass algorithm, as we have to read the hammer load file once to get the base value and the counts, and then
 * again to produce the modified file.
 *
 * The positional parameter is the name of the source hammer load file.  The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for enw hammer load file (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class RoleWeightProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleWeightProcessor.class);
    /** genome -> role -> count */
    private Map<String, CountMap<String>> countMaps;
    /** feature ID input column index */
    private int fidColIdx;
    /** role ID input column index */
    private int roleColIdx;
    /** mean group size */
    private double meanSize;

    // COMMAND-LINE OPTIONS

    @Argument(index = 0, metaVar = "hammers.tbl", usage = "input hammer load file", required = true)
    private File hammerFile;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure we can open the input file.
        if (! this.hammerFile.canRead())
            throw new FileNotFoundException("Input hammer file " + this.hammerFile + " is not found or unreadable.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Create the counting hash.
        this.countMaps = new HashMap<String, CountMap<String>>();
        // Open the input file for the counting pass.
        log.info("Processing {} for initial counts.", this.hammerFile);
        try (TabbedLineReader hammerStream = new TabbedLineReader(this.hammerFile)) {
            // Get the two critical input columns-- feature ID and role.
            this.fidColIdx = hammerStream.findField("fid");
            this.roleColIdx = hammerStream.findField("role");
            // Set up for progress messages.
            int inCount = 0;
            long lastMessage = System.currentTimeMillis();
            // Loop through the input file, updating counters.
            for (var line : hammerStream) {
                CountMap<String>.Count counter = this.getCounter(line);
                counter.increment();
                inCount++;
                if (log.isInfoEnabled()) {
                    long time = System.currentTimeMillis();
                    if (time - lastMessage >= 5000) {
                        log.info("{} hammers processed, {} genomes found.", inCount, this.countMaps.size());
                        lastMessage = time;
                    }
                }
            }
            log.info("Computing mean group size.");
            double total = 0.0;
            int nGroups = 0;
            for (var countMap : countMaps.values()) {
                for (var count : countMap.counts()) {
                    total += count.getCount();
                    nGroups++;
                }
            }
            if (nGroups == 0)
                throw new IOException("No hammers found in load file.");
            else
                this.meanSize = total / nGroups;
            log.info("Mean group size for {} groups is {}.", nGroups, this.meanSize);
        }
        log.info("Producing new hammer load file with updated weights.");
        try (TabbedLineReader hammerStream = new TabbedLineReader(this.hammerFile)) {
            // Get the weight column index.
            int weightColIdx = hammerStream.findField("strength");
            // Prime the output file with the headers.
            writer.println(hammerStream.header());
            // Set up for progress messages.
            int inCount = 0;
            long lastMessage = System.currentTimeMillis();
            // Loop through the input file again, producing output.
            for (var line : hammerStream) {
                inCount++;
                // Get the counter for this line's role/genome group.
                var counter = this.getCounter(line);
                // Compute the new weight.
                double weight = this.meanSize / counter.getCount();
                // Put it in the output.
                String[] fields = line.getFields();
                fields[weightColIdx] = Double.toString(weight);
                writer.println(StringUtils.join(fields, '\t'));
                if (log.isInfoEnabled()) {
                    long time = System.currentTimeMillis();
                    if (time - lastMessage >= 5000) {
                        log.info("{} hammer records copied.", inCount);
                        lastMessage = time;
                    }
                }
            }
        }
    }

    /**
     * @return the counter for the specified input line
     *
     * @param line	input line of interest
     *
     * @return the counter entry for the genome and role in the input line.
     */
    private CountMap<String>.Count getCounter(Line line) {
        String genomeId = Feature.genomeOf(line.get(this.fidColIdx));
        CountMap<String> map = this.countMaps.computeIfAbsent(genomeId, x -> new CountMap<String>());
        var retVal = map.getCounter(line.get(this.roleColIdx));
        return retVal;
    }

}
