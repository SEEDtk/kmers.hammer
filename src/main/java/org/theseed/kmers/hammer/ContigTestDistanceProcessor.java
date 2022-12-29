/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.distance.methods.AniDistanceMethod;
import org.theseed.genome.distance.methods.Measurer;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;
import org.theseed.utils.StringPair;

/**
 * This command analyzes the genome-level results from a hammer test against a synthetic sample (as produced by
 * "ContigTestProcessor").  The expected and actual genomes are computed, then the ANI distance calculated from
 * the test genome to both (or to the actual if they are the same).  We wish to determine (1) the distance range
 * to the actual genomes and (2) whether the actual genome is closer than the expected one.  This tells us how
 * good the hammers are at finding a close genome.
 *
 * The positional parameters are the name of the representative-genome directory (which contains all the genomes
 * that can be picked as expected or actual), the name of the source representative-genome database, and the name
 * of the test-genome directory.  Both of these should be GTO directories with standard naming (that is, the GTO
 * name is the genome ID with an extension of ".gto").
 *
 * The standard input should be the output of the hammer test, and the report will be generated on the standard output.
 *
 * The output contains the following flag columns.
 *
 * close		actual genome is within the maximum good-match distance
 * good			actual genome is the expected genome
 * better		actual genome is closer than the expected genome
 * near_miss	actual genome is only slightly further than the expected genome
 * alt_miss		actual genome is a qualified representative of the test genome
 * hard_miss	actual genome is much further than the expected genome
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing the test results (if not STDIN)
 * -o	output file to contain the report (if not STDOUT)
 *
 * --maxDist	maximum distance for an actual hit to be considered a good match (default 0.1)
 * --nearDist	maximum error for a miss to be considered a near-miss (default 0.01)
 * --cache		name of a file where genome ANI distances will be cached (default "distCache.tbl" in the current directory)
 * --method		method to use for determining actual hit (default STRENGTH)
 *
 * @author Bruce Parrello
 *
 */
public class ContigTestDistanceProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ContigTestDistanceProcessor.class);
    /** directory of test genomes */
    private GenomeDirectory testGenomes;
    /** directory of representative genomes */
    private GenomeDirectory repGenomes;
    /** ANI distance enging */
    private AniDistanceMethod distanceComputer;
    /** input column for test-genome hit location */
    private int locColIdx;
    /** input column for hitting hammer fid */
    private int hammerColIdx;
    /** input column for expected representative-genome ID */
    private int expectColIdx;
    /** input column for hammer strength */
    private int strengthColIdx;
    /** map of test genome IDs to score totals */
    private Map<String, WeightMap> scoreMap;
    /** map of test genome IDs to expected representatives */
    private Map<String, String> expectMap;
    /** cache of genome distances */
    private Map<StringPair, Double> distCache;
    /** distance measurer for the test genome */
    private Measurer testAnalysis;
    /** print writer for the distance cache file */
    private PrintWriter cacheStream;
    /** representative-genome database */
    private RepGenomeDb repDb;


    // COMMAND-LINE OPTIONS

    /** maximum distance for a closest-genome to be considered a good match */
    @Option(name = "--maxDist", aliases = { "-d" }, metaVar = "0.1", usage = "maximum distance for a match to be considered good")
    private double maxDist;

    /** maximum error for a closest-genome to be considered a near-miss */
    @Option(name = "--nearDist", metaVar = "0.05", usage = "maximum error for a match to be considered a near-miss")
    private double nearDist;

    /** flat file to use for cached distances */
    @Option(name = "--cache", metaVar = "distCache.tbl", usage = "flat file to contain cached genome distances")
    private File cacheFile;

    /** method for scoring hammer hits */
    @Option(name = "--method", usage = "method to use for scoring hammer hits")
    private HammerDb.Method scoreMethod;

    /** directory of representative genomes */
    @Argument(index = 0, metaVar = "repgenDir", usage = "directory of representative-genome GTOs",
            required = true)
    private File repDir;

    /** directory of test genomes */
    @Argument(index = 1, metaVar = "testgenDir", usage = "directory of test genome GTOs",
            required = true)
    private File testDir;

    /** original representative-genome database file */
    @Argument(index = 2, metaVar = "repDb.ser", usage = "representative-genome database file",
            required = true)
    private File repDbFile;

    @Override
    protected void setPipeDefaults() {
        this.maxDist = 0.05;
        this.nearDist = 0.01;
        this.cacheFile = new File(System.getProperty("user.dir"), "distCache.tbl");
        this.scoreMethod = HammerDb.Method.STRENGTH;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify the distance limits.
        if (this.maxDist < 0 || this.maxDist >= 1.0)
            throw new ParseFailureException("Maximum good-match distance must be between 0 and 1.");
        if (this.nearDist < 0 || this.nearDist >= 1.0)
            throw new ParseFailureException("Maximum near-miss distance must be between 0 and 1.");
        // Set up the rep-genome directory.
        if (! this.repDir.isDirectory())
            throw new FileNotFoundException("Representative-genome GTO directory " + this.repDir + " is not found or invalid.");
        this.repGenomes = new GenomeDirectory(this.repDir);
        // Load the rep-genome database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Representative-genome database file " + this.repDbFile + " is not found or unreadable.");
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // Set up the test-genome directory.
        if (! this.testDir.isDirectory())
            throw new FileNotFoundException("Test-genome GTO directory " + this.testDir + " is not found or invalid.");
        this.testGenomes = new GenomeDirectory(this.testDir);
        // Initialize the ANI-distance engine with default parameters.
        this.distanceComputer = new AniDistanceMethod();
        this.distanceComputer.parseParmString("");
        // Create the map from test-genome IDs to hit scores.
        final int hashSize = this.testGenomes.size() * 4 / 3 + 1;
        this.scoreMap = new HashMap<String, WeightMap>(hashSize);
        // Create the map for the expected genomes.
        this.expectMap = new HashMap<String, String>(hashSize);
        // Initialize the distance cache.
        this.distCache = new HashMap<StringPair, Double>(2000);
        if (! this.cacheFile.exists()) {
            // Here we must create the cache file.
            try (var writer = new PrintWriter(this.cacheFile)) {
                writer.println("genome1\tgenome2\tdistance");
            }
        } else {
            // Load the cache file into memory.
            try (TabbedLineReader cacheReader = new TabbedLineReader(this.cacheFile)) {
                for (var line : cacheReader) {
                    StringPair gPair = new StringPair(line.get(0), line.get(1));
                    double dist = line.getDouble(2);
                    this.distCache.put(gPair, dist);
                }
            }
            log.info("{} distance pairs loaded into cache from {}.", this.distCache.size(), this.cacheFile);
        }
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Get the input column indices.
        this.locColIdx = inputStream.findField("location");
        this.hammerColIdx = inputStream.findField("hammer_fid");
        this.expectColIdx = inputStream.findField("rep_id");
        this.strengthColIdx = inputStream.findField("strength");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Open the cache file for appending.
        this.cacheStream = null;
        try {
            var appender = new FileWriter(this.cacheFile, true);
            this.cacheStream = new PrintWriter(appender);
            // Get the list of test genomes.
            Set<String> testGenomeIDs = this.testGenomes.getGenomeIDs();
            // Create the count maps.  This insures each test genome has counts, even if it has no hits.
            for (String testGenomeId : testGenomeIDs)
                this.scoreMap.put(testGenomeId, new WeightMap());
            // First, we scan the input, counting hits.
            log.info("Scanning input file.");
            int inCount = 0;
            for (var line : inputStream) {
                String testGenomeId = StringUtils.substringBefore(line.get(this.locColIdx), ":");
                // Verify this line is valid.
                WeightMap counters = this.scoreMap.get(testGenomeId);
                if (counters == null)
                    throw new IOException("Invalid test-genome ID " + testGenomeId + " in input file.");
                String hitFid = line.get(this.hammerColIdx);
                String hitGenomeId = Feature.genomeOf(hitFid);
                String repId = line.get(this.expectColIdx);
                double strength = line.getDouble(this.strengthColIdx);
                // Insure the test genome is associated with its representative.
                this.expectMap.put(testGenomeId, repId);
                // Count the hit.
                HammerDb.Source hitSource = new HammerDb.Source(hitFid, strength);
                counters.count(hitGenomeId, this.scoreMethod.getWeight(hitSource));
                inCount++;
                if (log.isInfoEnabled() && inCount % 5000 == 0)
                    log.info("{} input lines processed.", inCount);
            }
            // Now, for each genome we look at the top genome (actual) and the expected genome and compute distances.
            // We need some counts.  This is the number of matches where expected == actual.
            int goodCount = 0;
            // This is the number of matches where the actual was within the minimum good-match distance.
            int closeCount = 0;
            // This is the number of matches where the actual was closer than the expected genome.
            int betterCount = 0;
            // This is the number of test genomes with no hits.
            int noMatchCount = 0;
            // This is the number of test genomes with a near-miss.
            int nearMatchCount = 0;
            // This is the number of test genomes matching an alternate representative.
            int altMatchCount = 0;
            // This is the number of test genomes with a bad match.
            int badMatchCount = 0;
            // Write the output header.
            writer.println("test_genome_id\ttest_genome_name\ttotal_hits\tmatch_genome_id\tmatch_hits\tmatch_dist\trep_id\trep_hits\trep_dist\tclose\tgood\tbetter\tnear_miss\talt_miss\thard_miss");
            // Loop through the test genomes.
            int gCount = 0;
            int gTotal = this.testGenomes.size();
            for (Genome testGenome : this.testGenomes) {
                gCount++;
                log.info("Processing genome {} of {}: {}.", gCount, gTotal, testGenome);
                String testGenomeId = testGenome.getId();
                String firstFields = testGenomeId + "\t" + testGenome.getName();
                // Denote we don't have a measurer yet.
                this.testAnalysis = null;
                // Get the best match.
                var counters = this.scoreMap.get(testGenomeId);
                var best = counters.getBestEntry();
                if (best == null) {
                    // Here there were no hits.  This almost NEVER happens.
                    writer.println(firstFields + "0\t\t\t\t\t\t\t\t\t\t\t\t");
                    noMatchCount++;
                } else {
                    // We will store the rating flags here.
                    String closeFlag = "";
                    String goodFlag = "";
                    String betterFlag = "";
                    String nearFlag = "";
                    String badFlag = "";
                    String altFlag = "";
                    // Get the genome IDs and hit counts.
                    String matchGenomeId = best.getKey();
                    double matchHits = best.getCount();
                    String repGenomeId = this.expectMap.get(testGenomeId);
                    double repHits = counters.getCount(repGenomeId);
                    double totalHits = counters.sum();
                    // Now compute the distances.
                    double matchDist = this.computeDistance(testGenome, matchGenomeId, "best");
                    if (matchDist <= this.maxDist) {
                        closeFlag = "Y";
                        closeCount++;
                    }
                    // Check to see if the expected genome is different.
                    double repDist;
                    if (repGenomeId.contentEquals(matchGenomeId)) {
                        // We matched the expected genome.
                        repDist = matchDist;
                        goodCount++;
                        goodFlag = "Y";
                    } else {
                        // Check the close representatives.
                        var closeList = this.repDb.findClose(testGenome);
                        boolean altRep = closeList.stream().anyMatch(x -> x.getGenomeId().contentEquals(matchGenomeId));
                        if (altRep) {
                            altFlag = "Y";
                            altMatchCount++;
                        }
                        // We need the distance to the expected genome.
                        repDist = this.computeDistance(testGenome, repGenomeId, "expected");
                        if (repDist >= matchDist) {
                            betterCount++;
                            betterFlag = "Y";
                        } else if (matchDist - repDist <= this.nearDist) {
                            nearMatchCount++;
                            nearFlag = "Y";
                        } else if (closeFlag.isEmpty() && altFlag.isEmpty()) {
                            // We matched the wrong genome and we are not close, so we count this as a hard miss.
                            badMatchCount++;
                            badFlag = "Y";
                        }
                    }
                    writer.println(firstFields + "\t" + totalHits + "\t" + matchGenomeId + "\t" + matchHits + "\t"
                            + matchDist + "\t" + repGenomeId + "\t" + repHits + "\t" + repDist + "\t" +
                            closeFlag + "\t" + goodFlag + "\t" + betterFlag + "\t" + nearFlag + "\t" +
                            altFlag + "\t" + badFlag);
                    writer.flush();
                }
            }
            // Write the final stats.
            log.info("{} test genomes processed.  {} good hits, {} close hits, {} better hits, {} with no hits, {} near misses, {} alternate representatives, {} hard misses.",
                    gTotal, goodCount, closeCount, betterCount, noMatchCount, nearMatchCount, altMatchCount, badMatchCount);
            log.info("Accuracy is {}.", (gCount - badMatchCount) / (double) gCount);
        } finally {
            // Note that closing the print writer automatically closes the underlying file writer.
            if (this.cacheStream != null)
                this.cacheStream.close();
        }
    }

    /**
     * Compute the distance between a test genome and the specified representative genome.
     *
     * @param testGenome	test genome from which the distance is computed
     * @param denomeId		ID of the other genome
     * @param string		display string for log messages
     *
     * @return the ANI distance between the two genomes
     */
    private double computeDistance(Genome testGenome, String genomeId, String string) {
        double retVal;
        StringPair gPair = new StringPair(testGenome.getId(), genomeId);
        Double dist = this.distCache.get(gPair);
        if (dist != null) {
            // Here we found the distance in the cache.
            retVal = dist;
        } else {
            // Here we must compute the distance.  We use the measurer and make sure the value is saved in the
            // cache file for the next time.
            if (this.testAnalysis == null)
                this.testAnalysis = this.distanceComputer.getMeasurer(testGenome);
            Genome otherGenome = this.repGenomes.getGenome(genomeId);
            log.info("Computing distance to {} match {}.", string, otherGenome);
            retVal = this.distanceComputer.getDistance(testAnalysis, otherGenome);
            this.cacheStream.println(gPair.getString1() + "\t" + gPair.getString2() + "\t" + Double.toString(retVal));
            this.cacheStream.flush();
        }
        return retVal;
    }

}
