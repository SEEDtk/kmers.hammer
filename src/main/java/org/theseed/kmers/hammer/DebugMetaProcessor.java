/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.slf4j.LoggerFactory;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerSourMap;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command analyzes the output from a synthetic sample contig test, filtering out the incorrect hits.  These
 * can be analyzed to determine how to improve results in generating classifiers.
 *
 * The standard input should contain the output of the contig test.  The filtered hits will be on the standard output.
 *
 * The positional parameters should be the name of the repXXX.stats.tbl file for the repgen set used to generate the
 * hammers, the name of the sour map save file for the hammers, and the name of the hammer load file.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing the contig-test output (if not STDIN)
 * -o	output file to contain the filtered hits (if not STDOUT)
 *
 * --counts		name of a file to contain a summary report (default "counts.tbl" in the current directory)
 *
 * @author Bruce Parrello
 *
 */
public class DebugMetaProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(DebugMetaProcessor.class);
    /** index of hammer-fid column */
    private int fidIdx;
    /** index of rep-id column */
    private int repIdIdx;
    /** index of strength column */
    private int strengthIdx;
    /** summary report data for each representative genome */
    private Map<String, RepInfo> repData;
    /** SOUR map for hammer fids */
    private HammerSourMap sourMap;
    /** summary report data for each hammer fid */
    private Map<String, FidInfo> fidData;


    // COMMAND-LINE OPTIONS

    /** output file for summary report */
    @Option(name = "--counts", metaVar = "summary.tbl", usage = "output file for summary report")
    private File countFile;

    /** name of the repXXX.stats.tbl file used to analyze the repgens found */
    @Argument(index = 0, metaVar = "repXX.stats.tbl", usage = "statistics file for the repgen set used to build the hammers",
            required = true)
    private File repStatsFile;

    /** name of the sour map save file for the hammer features */
    @Argument(index = 1, metaVar = "sourMap.ser", usage = "hammer SOUR map save file for the hammer features", required = true)
    private File sourMapFile;

    /** name of the hammer load file */
    @Argument(index = 2, metaVar = "hammers.tbl", usage = "hammer load file", required = true)
    private File hammerFile;

    /**
     * This object describes the basic data we need about each repgen genome.
     */
    private static class RepInfo {

        /** name of genome */
        private String name;
        /** number of neighbors */
        private int neighborhoodSize;

        /**
         * Construct a representative-info object from a rep-stats file input line.
         *
         * @param line		input data line from rep-stats file
         */
        protected RepInfo(TabbedLineReader.Line line) {
            this.name = line.get(1);
            this.neighborhoodSize = line.getInt(3);
        }

    }

    /**
     * This object tracks the information we need about each hammer fid.  It is sorted by the number
     * of bad hits (decreasing) and then the number of good hits (increasing).
     */
    private static class FidInfo implements Comparable<FidInfo> {

        /** hammer feature ID */
        private String fid;
        /** stats on strength for good hits */
        private SummaryStatistics goodHits;
        /** stats on strength for bad hits */
        private SummaryStatistics badHits;
        /** number of hammers for this fid */
        private int hammerCount;

        /**
         * Create a summary-report data holder for a hammer feature.
         *
         * @param fid	ID of the hammer feature of interest
         */
        protected FidInfo(String fid) {
            this.fid = fid;
            this.goodHits = new SummaryStatistics();
            this.badHits = new SummaryStatistics();
            this.hammerCount = 0;
        }

        @Override
        public int compareTo(FidInfo o) {
            int retVal = Long.compare(o.badHits.getN(), this.badHits.getN());
            if (retVal == 0) {
                retVal = Long.compare(this.goodHits.getN(), o.goodHits.getN());
                if (retVal == 0)
                    retVal = this.fid.compareTo(o.fid);
            }
            return retVal;
        }

    }

    @Override
    protected void setPipeDefaults() {
        this.countFile = new File(System.getProperty("user.dir"), "counts.tbl");
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Get the necessary column indices.  This also insures the input is compatible.
        this.fidIdx = inputStream.findField("hammer_fid");
        this.repIdIdx = inputStream.findField("rep_id");
        this.strengthIdx = inputStream.findField("strength");
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Read in the rep-stats file.
        this.repData = new HashMap<String, RepInfo>(4000);
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + " is not found or unreadable.");
        try (TabbedLineReader repStream = new TabbedLineReader(this.repStatsFile)) {
            if (repStream.findField("rep_name") != 1 || repStream.findField("members") != 3)
                throw new IOException("Repgen stats file " + this.repStatsFile + " is not the proper format.");
            int idColIdx = repStream.findField("rep_id");
            log.info("Reading repgen data from {}.", this.repStatsFile);
            for (var line : repStream) {
                String repId = line.get(idColIdx);
                this.repData.put(repId, new RepInfo(line));
            }
            log.info("{} repgen genomes found.", this.repData.size());
        }
        // Read in the sour map.
        if (! this.sourMapFile.canRead())
            throw new FileNotFoundException("Hammer sour map file " + this.sourMapFile + " is not found or unreadable.");
        this.sourMap = HammerSourMap.load(this.sourMapFile);
        // Create the fid map.
        this.fidData = new HashMap<String, FidInfo>(this.sourMap.size() * 4 / 3 + 1);
        // Validate the hammer load file.
        if (! this.hammerFile.canRead())
            throw new FileNotFoundException("Hammer load file " + this.hammerFile + " is not found or unreadable.");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Start by copying the header line.
        writer.println(inputStream.header());
        // Now loop through each input line.  We keep the ones where the genome ID for the hammer fid does NOT match the
        // expected rep-id.  We count the good and bad hits for each found genome ID.
        int keepCount = 0;
        int inCount = 0;
        long timer = System.currentTimeMillis();
        for (var line : inputStream) {
            inCount++;
            String fid = line.get(fidIdx);
            String foundRepId = Feature.genomeOf(fid);
            String expectedRepId = line.get(repIdIdx);
            RepInfo foundInfo = this.repData.get(foundRepId);
            if (foundInfo == null)
                throw new IOException("Invalid rep ID " + foundRepId + " found in input file.");
            double strength = line.getDouble(strengthIdx);
            FidInfo fidInfo = this.fidData.computeIfAbsent(fid, x -> new FidInfo(x));
            if (! foundRepId.contentEquals(expectedRepId)) {
                keepCount++;
                writer.println(line.toString());
                fidInfo.badHits.addValue(strength);
            } else
                fidInfo.goodHits.addValue(strength);
            if (log.isInfoEnabled() && System.currentTimeMillis() - timer >= 5000) {
                timer = System.currentTimeMillis();
                log.info("{} lines read.  {} kept.", inCount, keepCount);
            }
        }
        log.info("{} lines read; {} kept.", inCount, keepCount);
        // Now get the hammer counts from the hammer load file.
        try (TabbedLineReader hammerStream = new TabbedLineReader(this.hammerFile)) {
            int hFidIdx = hammerStream.findField("fid");
            log.info("Scanning hammer file {} to count hammers per fid.", this.hammerFile);
            for (var line : hammerStream) {
                String fid = line.get(hFidIdx);
                FidInfo fidInfo = this.fidData.computeIfAbsent(fid, x -> new FidInfo(x));
                fidInfo.hammerCount++;
            }
        }
        // Write out the summary report.
        var counts = new ArrayList<FidInfo>(this.fidData.values());
        Collections.sort(counts);
        try (PrintWriter summaryWriter = new PrintWriter(this.countFile)) {
            summaryWriter.println("hammer_fid\trole\trep_name\tneighbor_size\thammer_count\tgood_hits\tbad_hits\tgood_frac\tmean_good\tmean_bad\tmin_good_strength\tmax_bad_strength");
            for (var count : counts) {
                long good = count.goodHits.getN();
                long bad = count.badHits.getN();
                if (good > 0 || bad > 0) {
                    var frac = good / (double) (good + bad);
                    String repId = this.sourMap.getRepId(count.fid);
                    String roleId = this.sourMap.getRole(count.fid);
                    RepInfo repInfo = this.repData.get(repId);
                    summaryWriter.println(count.fid + "\t" + roleId + "\t" + repInfo.name + "\t" + repInfo.neighborhoodSize + "\t"
                            + count.hammerCount + "\t" + + good + "\t" + bad + "\t" + frac + "\t" + count.goodHits.getMean() + "\t"
                            + count.badHits.getMean() + "\t" + count.goodHits.getMin() + "\t" + count.badHits.getMax());
                }
            }
        }
    }

}
