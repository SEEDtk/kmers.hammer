/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command looks at the results from a contig test performed by ContigTestProcessor on a synthetic sample
 * produced by "SyntheticSampleProcessor" (in "genome.bins") and determines how well the hammers did in determining
 * the correct genome on a contig and genome basis.
 *
 * The input file is tab-delimited with five columns:  "location" is the contig location of the hammer hit, "hammer_fid"
 * is the feature ID of the hammer that hit the contig, "genome_name" is the name of the genome that contained the contig,
 * "rep_id" is the expected representative genome, and "distance" is the distance from the genome that contained the contig
 * to that representative.  The contig ID in the location will be prefixed by the genome ID with a colon (:).  The input
 * file is sorted by location, so therefore also by contig within genome.
 *
 * For each contig, we determine the representative genomes that hit it and how many hits there were, along with the ID of
 * the representative that hit the most.  We also summarize this information for each whole genome as well.
 *
 * The positional parameter is the name of a repXXX.stats.tbl file containing the genome names for the representing genomes.
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	name of the input file containing the test results (if not STDIN)
 * -o	name of the output file for the report (if not STDOUT)
 *
 * --summary	if specified, only genome records will be output
 * --limit		maximum number of hit targets to output per genome (default 10)
 *
 * @author Bruce Parrello
 *
 */
public class ContigTestAnalysisProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ContigTestAnalysisProcessor.class);
    /** index (0-based) of input column containing hit location */
    private int locationColIdx;
    /** index (0-based) of input column containing hammer ID */
    private int hitColIdx;
    /** index (0-based) of input column containing contig genome name */
    private int nameColIdx;
    /** index (0-based) of input column containing representing genome ID */
    private  int repColIdx;
    /** index (0-based) of input column containing distance */
    private int distIdx;
    /** map of representing genome IDs to names */
    private Map<String, String> rNameMap;

    // COMMAND-LINE OPTIONS

    /** TRUE to disable detail reports */
    @Option(name = "--summary", usage = "if specified, only genome-level statistics will be output")
    private boolean summaryOnly;

    /** maxmimum number of hit counts to output per object */
    @Option(name = "--limit", usage = "maximum number of hit counts to output per object (0 for all)")
    private int limitCount;

    /** stats file containing representative genome IDs and names */
    @Argument(index = 0, metaVar = "repXX.stats.tbl", usage = "name of a file containing representing genome IDs and names")
    private File repFileName;

    /**
     * This class contains useful data about each genome.  it is sorted by expected-hit ratio, meaning the worst genomes
     * float to the top.
     */
    private static class GenomeData implements Comparable<GenomeData> {

        /** ID of the genome */
        private String genomeId;
        /** name of the genome */
        private String genomeName;
        /** ID of the expected representative */
        private String repId;
        /** distance to the expected representative */
        private double distance;
        /** map of contig IDs to hit counts for each repId */
        private Map<String, CountMap<String>> contigMap;
        /** hit counts for each rep ID for the whole genome */
        private CountMap<String> counts;
        /** number of total hits */
        private int totalHits;
        /** ratio of expected hits to total hits */
        private double ratio;

        /**
         * Create a new genome-data object for the specified genome.
         *
         * @param id		ID of originating genome
         * @param name		name of originating genome
         * @param rep_id	ID of expected representative
         * @param dist		distance to expected representative
         */
        protected GenomeData(String id, String name, String rep_id, double dist) {
            // Save the basic data.
            this.genomeId = id;
            this.genomeName = name;
            this.repId = rep_id;
            this.distance = dist;
            // Initialize the counters.
            this.counts = new CountMap<String>();
            this.contigMap = new TreeMap<String, CountMap<String>>();
            this.totalHits = 0;
            this.ratio = 0.0;
        }

        /**
         * Count the specified representative in the specified contig.
         *
         * @param contigID	ID of the contig
         * @param hit_rep	ID of the representative found
         */
        protected void count(String contigID, String hit_rep) {
            this.counts.count(hit_rep);
            CountMap<String> counter = this.contigMap.computeIfAbsent(contigID, x -> new CountMap<String>());
            counter.count(hit_rep);
            // Update the ratio.
            this.totalHits++;
            double num = this.counts.getCount(this.repId);
            this.ratio = num / this.totalHits;
        }

        /**
         * @return the originating genome name
         */
        protected String getGenomeName() {
            return this.genomeName;
        }

        /**
         * @return the ID of the expected representative
         */
        protected String getRepId() {
            return this.repId;
        }

        /**
         * @return the distance to the expected representative
         */
        protected double getDistance() {
            return this.distance;
        }

        /**
         * @return the map of contig IDs to contig hit counters
         */
        protected Map<String, CountMap<String>> getContigMap() {
            return this.contigMap;
        }

        /**
         * @return the genome hit counter
         */
        protected CountMap<String> getCounts() {
            return this.counts;
        }

        @Override
        public int compareTo(GenomeData o) {
            int retVal = Double.compare(this.ratio, o.ratio);
            if (retVal == 0) {
                // If two ratios are equal, the genome closest to its representative is worse.
                retVal = Double.compare(this.distance, o.distance);
                // In an emergency we sort by genome name.
                if (retVal == 0)
                    retVal = this.genomeName.compareTo(o.genomeName);
            }
            return retVal;
        }

        /**
         * @return the ID of the originating genome
         */
        public String getGenomeID() {
            return this.genomeId;
        }

    }

    @Override
    protected void setPipeDefaults() {
        this.summaryOnly = false;
        this.limitCount = 10;
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Find the input file fields.
        this.locationColIdx = inputStream.findField("location");
        this.hitColIdx = inputStream.findField("hammer_fid");
        this.repColIdx = inputStream.findField("rep_id");
        this.nameColIdx = inputStream.findField("genome_name");
        this.distIdx = inputStream.findField("distance");
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify the stats file name.
        if (! this.repFileName.canRead())
            throw new FileNotFoundException("Stats file name " + this.repFileName + " is not found or unreadable.");
        // Read in the IDs and names.
        this.rNameMap = TabbedLineReader.readMap(this.repFileName, "rep_id", "rep_name");
        log.info("{} representing genome IDs found in {}.", this.rNameMap.size(), this.repFileName);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // This maps each incoming genome ID to the useful data we need about it.
        Map<String, GenomeData> gDataMap = new HashMap<String, GenomeData>(1000);
        // Loop through the input, counting each hit.
        log.info("Reading hits from input.");
        int hitCount = 0;
        for (var line : inputStream) {
            Location loc = Location.parseSeedLocation(line.get(this.locationColIdx));
            String contigId = loc.getContigId();
            String genomeId = StringUtils.substringBefore(contigId, ":");
            String repId = line.get(this.repColIdx);
            double dist = line.getDouble(this.distIdx);
            String name = line.get(this.nameColIdx);
            // Get the descriptor for the genome hit.
            GenomeData gData = gDataMap.computeIfAbsent(genomeId, x -> new GenomeData(x, name, repId, dist));
            // Compute the genome hit.
            String hitGenomeId = Feature.genomeOf(line.get(this.hitColIdx));
            // Record the hit.
            gData.count(contigId, hitGenomeId);
            hitCount++;
            if (log.isInfoEnabled() && hitCount % 100000 == 0)
                log.info("{} hits processed.", hitCount);
        }
        log.info("{} hits processed from {} genomes.", hitCount, gDataMap.size());
        // Now we write the output.
        log.info("Producing report.");
        writer.println("object_id\tgenome_name\thit_id\thit_name\tcount\texpected\tdist_to_exp");
        // Loop through the genome descriptors in natural sort order.
        gDataMap.values().stream().sorted().forEach(x -> this.writeGenome(writer, x));
    }

    /**
     * Write the data for a specified genome.
     *
     * @param writer		output writer for the report
     * @param gDatum		data collected on the originating genome
     */
    private void writeGenome(PrintWriter writer, GenomeData gDatum) {
        // Output a spacer for readability.
        writer.println();
        // Start with the genome summary.
        var counts = gDatum.getCounts();
        this.showHitCounts(writer, gDatum.getGenomeID(), gDatum, counts);
        if (! this.summaryOnly) {
            // Now write the hit counts for each contig.
            var contigMap = gDatum.getContigMap();
            contigMap.entrySet().stream().forEach(x -> this.showHitCounts(writer, x.getKey(), gDatum, x.getValue()));
        }
        log.info("{} ({}) written to output.", gDatum.getGenomeID(), gDatum.getGenomeName());
    }

    /**
     * Write out the hit counts for a genome or contig.
     * @param writer		output print writer
     * @param sourceId		ID of object relevant to the counts
     * @param gDatum		descriptor of the object's genome
     * @param counts		map of representing genome IDs to hit counts
     */
    private void showHitCounts(PrintWriter writer, String sourceId, GenomeData gDatum, CountMap<String> counts) {
        // Loop through the hit counts.  Note we stop at the limit count.
        List<CountMap<String>.Count> sortedCounts = counts.sortedCounts();
        int n = sortedCounts.size();
        if (this.limitCount > 0 && this.limitCount < n)
            n = this.limitCount;
        var countView = sortedCounts.subList(0, n);
        for (var count : countView) {
            String targetId = count.getKey();
            String targetName = this.rNameMap.getOrDefault(targetId, "");
            // Output this hit count.
            writer.format("%s\t%s\t%s\t%s\t%d\t%s\t%6.2f%n", sourceId, gDatum.getGenomeName(), targetId, targetName,
                    count.getCount(), gDatum.getRepId(), gDatum.getDistance());
        }
    }


}
