/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.utils.BasePipeProcessor;

/**
 * This command analyzes the results from a hammer contig test to determine useful statistics.  For each test-genome/rep-genome
 * pairing for which at least one hammers hit, it will output the number of hits, the total strength, the size of the
 * rep-genome neighborhood used during hammer computation, and whether or not the match was acceptable.  To do this, it
 * needs the repgen list file used to create the hammers, the output of the ContigTestProcessor (which lists every single hit),
 * and the output of ContigTestDistanceProcessor (which tells us which hits are good).  All this information is collated
 * to produce a single report.
 *
 * The positional parameters are the name of the repgen list file (usually of the form "repXXX.list.tbl"), the name of the
 * actual repgen database file (usually of the form "repXXX.ser"), the name of a genome source for the test genomes, and the
 * name of a file that tells us which matches are good.  This latter file is tab-delimited, and we need the test-genome ID
 * and name (columns "test_genome_id" and "test_genome_name"), the actual match ID ("match_genome_id"), a flag telling
 * us if the actual match is bad ("hard_miss"), and a flag telling us if the strongest match is good ("strength_good"). This
 * conforms to the output of ContigTestDistanceProcessor.
 *
 * The standard input should contain the output of ContigTestProcessor, which contains details on every single hammer hit from
 * the test.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	name of the input file (if not STDIN)
 * -o	name of the output file for the report (if not STDOUT)
 *
 * --type		type of source for the test genomes (default DIR)
 * --summary	if specified, only hit counts for actual or acceptable matches will be output
 *
 * @author Bruce Parrello
 *
 */
public class ContigTestStatisticsProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ContigTestStatisticsProcessor.class);
    /** map of reference-genome IDs to neighbor counts */
    private CountMap<String> neighborhoodCounts;
    /** map of test genome IDs to match statistics */
    private Map<String, TestGenomeData> testGenomeMap;
    /** index of the hit-location column in the input file */
    private int locColIdx;
    /** index of the hammer-id location column in the input file */
    private int hammerColIdx;
    /** index of the strength column in the input file */
    private int strengthColIdx;

    // COMMAND-LINE OPTIONS

    /** test genome source type */
    @Option(name = "--type", usage = "test genome source type")
    private GenomeSource.Type testSourceType;

    /** if specified, only summary lines will be output */
    @Option(name = "--summary", usage = "if specified, only lines with actual or acceptable matches will be output")
    private boolean summaryFlag;

    /** representative-genome list file containing neighborhood information */
    @Argument(index = 0, metaVar = "repXX.list.tbl", usage = "name of the repgen set list file used to generate the hammers",
            required = true)
    private File repListFile;

    /** representative-genome database file used to determine good matches */
    @Argument(index = 1, metaVar = "repXX.ser", usage = "name of the repgen database for the repgen set used to generate the hammers",
            required = true)
    private File repDbFile;

    /** test genome source file or directory */
    @Argument(index = 2, metaVar = "testGenomeDir", usage = "genome source for the test genomes")
    private File testGenomeDir;

    /** distance-analysis file showing hard misses */
    @Argument(index = 3, metaVar = "contigTest.dists.tbl", usage = "contig-test distance analysis file")
    private File distFile;

    /**
     * This utility class describes all the data we need about a match between
     * a test genome and a representative.  We sort from most hits to least.
     */
    protected static class TestMatchData implements Comparable<TestMatchData> {

        /** genome ID of the matched genome */
        private String genomeId;
        /** TRUE if the match is acceptable */
        private boolean acceptable;
        /** number of hammer hits */
        private int hits;
        /** total strength of all hammer hits */
        private double strengthTotal;

        /**
         * Create a new, blank match descriptor.
         *
         * @param id			ID of the matched genome
         */
        public TestMatchData(String id) {
            this.genomeId = id;
            this.acceptable = false;
            this.hits = 0;
            this.strengthTotal = 0.0;
        }

        @Override
        public int compareTo(TestMatchData o) {
            int retVal = o.hits - this.hits;
            // If two hit counts are the same, pick the one with the highest worth.
            if (retVal == 0) {
                retVal = Double.compare(o.strengthTotal, this.strengthTotal);
                // If both counters are the same, sort by match genome ID.
                if (retVal == 0)
                    retVal = this.genomeId.compareTo(o.genomeId);
            }
            return retVal;
        }

    }

    /**
     * This utility class describes all the data we need about a test genome.
     */
    protected static class TestGenomeData implements Comparable<TestGenomeData>, Iterable<TestMatchData> {

        /** genome ID */
        private String genomeId;
        /** genome name */
        private String genomeName;
        /** ID of actual match */
        private String actualId;
        /** ID of strongest match */
        private String strongId;
        /** map of match genome IDs to match statistics */
        private Map<String, TestMatchData> matchMap;
        /** error margin between the actual match and the first acceptable match */
        private int errorMargin;
        /** sorted list of match pairing objects */
        private List<TestMatchData> matchList;

        /**
         * Create a blank test-genome descriptor.
         *
         * @param id		ID of the test genome
         * @param name		name of the test genome
         */
        public TestGenomeData(String id, String name) {
            this.genomeId = id;
            this.genomeName = name;
            this.matchMap = new HashMap<String, TestMatchData>(20);
            this.actualId = "";
            this.strongId = "";
        }

        /**
         * Record a possible match as bad or acceptable
         *
         * @param matchId		ID of the matching genome
         * @param goodFlag		TRUE if the match is acceptable
         */
        public void addMatch(String matchId, boolean goodFlag) {
            TestMatchData matchData = findMatch(matchId);
            matchData.acceptable = goodFlag;
        }

        /**
         * Find the match object for a given matching genome ID.  If it does not exist it will be created.
         *
         * @param matchId	the ID of the genome of interest
         *
         * @return the match object for the identified genome
         */
        public TestMatchData findMatch(String matchId) {
            TestMatchData matchData = this.matchMap.computeIfAbsent(matchId, x -> new TestMatchData(x));
            return matchData;
        }

        /**
         * Denote that a particular match is the actual match.
         *
         * @param matchId	ID of the matched genome
         */
        public void setActual(String matchId) {
            // Save the actual-match ID.
            this.actualId = matchId;
            // Insure the pairing is set up.
            this.findMatch(matchId);
        }

        /**
         * Denote that a particular match is the strongest match.
         *
         * @param matchId	ID of the matched genome
         */
        public void setStrongest(String matchId) {
            // Save the actual-match ID.
            this.strongId = matchId;
            // Insure the pairing is set up.
            this.findMatch(matchId);
        }

        /**
         * Record a hammer hit against the specified match genome.
         *
         * @param matchGenomeId		representative genome for the hit
         * @param strength			strength of the hammer
         */
        public void recordHit(String matchGenomeId, double strength) {
            var matchData = this.findMatch(matchGenomeId);
            matchData.hits++;
            matchData.strengthTotal += strength;
        }

        /**
         * Compute the error margin and sort the pairings.  This should be
         * done after all the hits are counted.
         */
        public void sortMatches() {
            this.matchList = this.matchMap.values().stream().sorted().collect(Collectors.toList());
            // The closest representative is always in the list, even if it got 0 hits.  So, there
            // is always an acceptable genome.  If the first genome is acceptable, the error margin
            // is 0.
            if (this.matchList.get(0).acceptable)
                this.errorMargin = 0;
            else {
                // Otherwise, we find the first acceptable genome and subtract its hit count from the
                // one with the most hits.
                int bestHits = this.matchList.get(0).hits;
                int goodHits = this.matchList.stream().filter(x -> x.acceptable).findFirst().get().hits;
                this.errorMargin = bestHits - goodHits;
            }
        }

        @Override
        public int compareTo(TestGenomeData o) {
            // Sort from highest error margin to lowest, then sort by genome ID.
            int retVal = o.errorMargin - this.errorMargin;
            if (retVal == 0)
                retVal = this.genomeId.compareTo(o.genomeId);
            return retVal;
        }

        @Override
        public Iterator<TestMatchData> iterator() {
            return this.matchList.iterator();
        }

    }

    @Override
    protected void setPipeDefaults() {
        this.testSourceType = GenomeSource.Type.DIR;
        this.summaryFlag = false;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Insure that the input files all exist.
        if (! this.distFile.canRead())
            throw new FileNotFoundException("Distance-analysis file " + this.distFile + " is not found or unreadable.");
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Repgen database file " + this.repDbFile + " is not found or unreadable.");
        if (! this.repListFile.canRead())
            throw new FileNotFoundException("Repgen list file " + this.repListFile + " is not found or unreadable.");
        if (! this.testGenomeDir.exists())
            throw new FileNotFoundException("Test genome source " + this.testGenomeDir + " is not found.");
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Verify that we have a real contig-test output file by retrieving the column indices.
        this.locColIdx = inputStream.findField("location");
        this.hammerColIdx = inputStream.findField("hammer_fid");
        this.strengthColIdx = inputStream.findField("strength");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // First we compute the hard misses and the good actual matches.
        this.getHardMisses();
        // Next we load the repgen database and process the test-genome source.  This gives us the names and IDs of the
        // test genomes, as well as adding to our list of acceptable matches.
        this.processRepGen();
        // Now we compute the neighborhood counts from the repgen list file.  We need these for our output.
        this.processRepList();
        // Finally, we are ready to count hits.
        log.info("Processing hammer hits from input.");
        int inCount = 0;
        for (var line : inputStream) {
            // Get the test genome ID from the location.
            String testGenomeId = StringUtils.substringBefore(line.get(this.locColIdx), ":");
            // Get the match genome ID from the hammer fid.
            String matchGenomeId = Feature.genomeOf(line.get(this.hammerColIdx));
            // Get the strength.
            double strength = line.getDouble(this.strengthColIdx);
            // Find the test genome.
            TestGenomeData testGenomeData = this.testGenomeMap.get(testGenomeId);
            if (testGenomeData == null)
                throw new IOException("Test genome " + testGenomeId + " not found in support files:  input may not be from the same run.");
            // Record the hit.
            testGenomeData.recordHit(matchGenomeId, strength);
            inCount++;
            if (log.isInfoEnabled() && inCount % 10000 == 0)
                log.info("{} input lines processed.", inCount);
        }
        log.info("{} hits counted.  Computing error margins.", inCount);
        this.testGenomeMap.values().stream().forEach(x -> x.sortMatches());
        // Sort the test genome data.
        log.info("Sorting {} count sets.", this.testGenomeMap.size());
        List<TestGenomeData> testCounts = this.testGenomeMap.values().stream().sorted().collect(Collectors.toList());
        // Now we write out all the hit totals.  We write an output line for each matched pair.
        log.info("Writing report.");
        writer.println("genome_id\tgenome_name\trep_id\thits\tstrength_total\tstrength_mean\tneighborhood\tacceptable\tactual\tstrongest");
        // Loop through the test genomes.
        for (var testGenomeData : testCounts) {
            log.info("Writing counts for {}: {}.", testGenomeData.genomeId, testGenomeData.genomeName);
            // Remember the test genome ID and name.
            String genomeId = testGenomeData.genomeId;
            String genomeName = testGenomeData.genomeName;
            // Remember the actual match ID and the strongest-match ID.
            String actualId = testGenomeData.actualId;
            String strongId = testGenomeData.strongId;
            // Loop through the match pairs.
            for (var testMatchData : testGenomeData) {
                // Determine if this is the actual matched genome or the strongest-match genome.
                String actualFlag = (testMatchData.genomeId.equals(actualId) ? "*" : "");
                String strongFlag = (testMatchData.genomeId.equals(strongId) ? "Y" : "");
                // Filter detail lines if this is a summary report.
                if (! this.summaryFlag || ! actualFlag.isEmpty() || testMatchData.acceptable
                        || ! strongFlag.isEmpty()) {
                    // Determine if the match is acceptable.
                    String goodFlag = (testMatchData.acceptable ? "Y" : "");
                    // Compute the mean worth.
                    double meanStrength;
                    if (testMatchData.hits <= 0)
                        meanStrength = 0.0;
                    else
                        meanStrength = testMatchData.strengthTotal / testMatchData.hits;
                    // Get the matched genome's neighborhood size.
                    int neighborCount = this.neighborhoodCounts.getCount(testMatchData.genomeId);
                    writer.println(genomeId + "\t" + genomeName + "\t" + testMatchData.genomeId + "\t"
                                + Integer.toString(testMatchData.hits) + "\t"
                                + Double.toString(testMatchData.strengthTotal) + "\t" +  Double.toString(meanStrength)
                                + "\t" + Integer.toString(neighborCount) + "\t" + goodFlag + "\t" + actualFlag
                                + "\t" + strongFlag);
                }
            }
        }
    }

    /**
     * This method reads the distance analysis file.  It creates the main test-genome map and gets us an initial list of
     * good matches and hard misses.
     *
     * @throws IOException
     */
    private void getHardMisses() throws IOException {
        // Create the master test-genome hash.
        this.testGenomeMap = new HashMap<String, TestGenomeData>(2000);
        // Process the distance analysis file.  For each test genome, we need to get its actual match and whether it was
        // good or bad.
        log.info("Scanning distance-analysis file {}.", this.distFile);
        try (TabbedLineReader inStream = new TabbedLineReader(this.distFile)) {
            // Locate the input columns.
            int genomeCol = inStream.findField("test_genome_id");
            int nameCol = inStream.findField("test_genome_name");
            int matchCol = inStream.findField("match_genome_id");
            int repCol = inStream.findField("rep_id");
            int missCol = inStream.findField("hard_miss");
            int strongCol = inStream.findField("strong_genome_id");
            int strengthFlagCol = inStream.findField("strength_good");
            // Set up some counters.
            int inCount = 0;
            int goodCount = 0;
            int badCount = 0;
            // Loop through the file.
            for (var line : inStream) {
                String genomeId = line.get(genomeCol);
                String matchId = line.get(matchCol);
                boolean hardMissFlag = line.getFlag(missCol);
                // Create the test-genome counter.
                TestGenomeData genomeData = new TestGenomeData(genomeId, line.get(nameCol));
                // Add the match data.  Note that a good match is whatever is NOT a hard miss,
                // and that this is the so-called "actual" match.
                genomeData.addMatch(matchId, ! hardMissFlag);
                genomeData.setActual(matchId);
                // Add the strongest-match data.  It is a good match if the strength-good
                // flag is set.
                String strongId = line.get(strongCol);
                boolean strongFlag = line.getFlag(strengthFlagCol);
                genomeData.addMatch(strongId, strongFlag);
                genomeData.setStrongest(strongId);
                // Update the counters.
                inCount++;
                if (hardMissFlag) {
                    badCount++;
                    // If we have a hard miss, add the expected value as a good match.
                    genomeData.addMatch(line.get(repCol), true);
                } else
                    goodCount++;
                // Put the test-genome counter in the map.
                this.testGenomeMap.put(genomeId, genomeData);
            }
            log.info("{} test genomes processed, {} good matches, {} hard misses.", inCount, goodCount, badCount);
        }
    }

    /**
     * This method processes the test genomes against the repgen database to get a comprehensive list of good matches
     * based on seed-protein similarity.  This nets us not only the expected matches, but matches for representatives
     * whose neighborhoods overlap the expected match's neighborhood.  Such matches are also acceptable.
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    private void processRepGen() throws IOException, ParseFailureException {
        log.info("Loading repgen database from {}.", this.repDbFile);
        RepGenomeDb repdb = RepGenomeDb.load(this.repDbFile);
        log.info("Loading test genomes from {}.", this.testGenomeDir);
        GenomeSource testGenomes = this.testSourceType.create(this.testGenomeDir);
        // This will be slow.  For each test genome, we find all the eligible representatives.
        log.info("Scanning for representation.");
        for (Genome testGenome : testGenomes) {
            Collection<RepGenomeDb.Representation> reps = repdb.findClose(testGenome);
            log.info("{} representatives found for {}.", reps.size(), testGenome);
            // We must mark the representatives as good matches for this test genome.
            TestGenomeData testGenomeData = this.testGenomeMap.computeIfAbsent(testGenome.getId(),
                    x -> new TestGenomeData(x, testGenome.getName()));
            for (var rep : reps) {
                testGenomeData.addMatch(rep.getGenomeId(), true);
            }
        }
        if (log.isInfoEnabled()) {
            // Now we run through all the matches, counting good and bad ones.
            int goodCount = 0;
            int badCount = 0;
            int possibleCount = 0;
            int gCount = 0;
            for (var testGenomeData : this.testGenomeMap.values()) {
                gCount++;
                String actualId = testGenomeData.actualId;
                for (var testGenomeMatch : testGenomeData.matchMap.values()) {
                    if (testGenomeMatch.acceptable) {
                        possibleCount++;
                        if (testGenomeMatch.genomeId.equals(actualId))
                            goodCount++;
                    } else if (testGenomeMatch.genomeId.equals(actualId))
                        badCount++;
                }
            }
            log.info("{} test genomes.  {} good matches, {} bad matches, {} acceptable match pairs.",
                    gCount, goodCount, badCount, possibleCount);
        }
    }

    /**
     * Read the representative-genome list file and determine the size of each representative genome's
     * neighborhood.
     *
     * @throws IOException
     */
    private void processRepList() throws IOException {
        log.info("Analyzing neighborhood data from {}.", this.repListFile);
        // Our basic strategy is to count the number of times each representative occurs in the "rep_id"
        // column.
        this.neighborhoodCounts = new CountMap<String>();
        try (TabbedLineReader inStream = new TabbedLineReader(this.repListFile)) {
            int repIdColIdx = inStream.findField("rep_id");
            int inCount = 0;
            for (var line : inStream) {
                String repId = line.get(repIdColIdx);
                this.neighborhoodCounts.count(repId);
                inCount++;
                if (log.isInfoEnabled() && inCount % 10000 == 0)
                    log.info("{} list-file lines processed.", inCount);
            }
            log.info("{} lines read from {}.  {} neighborhoods found.", inCount, this.repListFile,
                    this.neighborhoodCounts.size());
        }
    }

}
