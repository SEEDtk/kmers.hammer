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

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.genome.distance.methods.AniDistanceMethod;
import org.theseed.genome.distance.methods.Measurer;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command analyzes the genome-level results from a hammer test against a synthetic sample (as produced by
 * "ContigTestProcessor").  The expected and actual genomes are computed, then the ANI distance calculated from
 * the test genome to both (or to the actual if they are the same).  We wish two determine (1) the distance range
 * to the actual genomes and (2) whether the actual genome is closer than the expected one.  This tells us how
 * good the hammers are at finding a close genome.
 *
 * The positional parameters are the name of the representative-genome directory (which contains all the genomes
 * that can be picked as expected or actual) and the name of the test-genome directory.  Both of these should be
 * GTO directories with standard naming (that is, the GTO name is the genome ID with an extension of ".gto").
 *
 * The standard input should be the output of the hammer test, and the report will be generated on the standard output.
 *
 * The output contains the following flag columns.
 *
 * close		actual genome is within the maximum good-match distance
 * good			actual genome is the expected genome
 * better		actual genome is closer than the expected genome
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing the test results (if not STDIN)
 * -o	output file to contain the report (if not STDOUT)
 *
 * --maxDist	maximum distance for an actual hit to be considered a good match
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
    /** map of test genome IDs to hit counts */
    private Map<String, CountMap<String>> hitCountMap;
    /** map of test genome IDs to expected representatives */
    private Map<String, String> expectMap;

    // COMMAND-LINE OPTIONS

    /** maximum distance for a closest-genome to be considered a good match */
    @Option(name = "--maxDist", aliases = { "-d" }, metaVar = "0.1", usage = "maximum distance for a match to be considered good")
    private double maxDist;

    /** directory of representative genomes */
    @Argument(index = 0, metaVar = "repgenDir", usage = "directory of representative-genome GTOs")
    private File repDir;

    /** directory of test genomes */
    @Argument(index = 1, metaVar = "testgenDir", usage = "directory of test genome GTOs")
    private File testDir;

    @Override
    protected void setPipeDefaults() {
        this.maxDist = 0.05;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Verify the distance limit.
        if (this.maxDist < 0 || this.maxDist >= 1.0)
            throw new ParseFailureException("Maximum good-match distance must be between 0 and 1.");
        // Set up the rep-genome directory.
        if (! this.repDir.isDirectory())
            throw new FileNotFoundException("Representative-genome GTO directory " + this.repDir + " is not found or invalid.");
        this.repGenomes = new GenomeDirectory(this.repDir);
        // Set up the test-genome directory.
        if (! this.testDir.isDirectory())
            throw new FileNotFoundException("Test-genome GTO directory " + this.testDir + " is not found or invalid.");
        this.testGenomes = new GenomeDirectory(this.testDir);
        // Initialize the ANI-distance engine with default parameters.
        this.distanceComputer = new AniDistanceMethod();
        this.distanceComputer.parseParmString("");
        // Create the map from test-genome IDs to hit counts.
        this.hitCountMap = new HashMap<String, CountMap<String>>(this.testGenomes.size() * 4 / 3 + 1);
        // Create the map for the expected genomes.
        this.expectMap = new HashMap<String, String>(this.testGenomes.size() * 4 / 3 + 1);
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Get the input column indices.
        this.locColIdx = inputStream.findField("location");
        this.hammerColIdx = inputStream.findField("hammer_fid");
        this.expectColIdx = inputStream.findField("rep_id");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Get the list of test genomes.
        Set<String> testGenomeIDs = this.testGenomes.getGenomeIDs();
        // Create the hit-count map.  This insures each test genome has counts, even if it has no hits.
        for (String testGenomeId : testGenomeIDs)
            this.hitCountMap.put(testGenomeId, new CountMap<String>());
        // First, we scan the input, counting hits.
        log.info("Scanning input file.");
        int inCount = 0;
        for (var line : inputStream) {
            String testGenome = StringUtils.substringBefore(line.get(this.locColIdx), ":");
            // Verify this line is valid.
            CountMap<String> counters = this.hitCountMap.get(testGenome);
            if (counters == null)
                throw new IOException("Invalid test-genome ID " + testGenome + " in input file.");
            String hitGenome = Feature.genomeOf(line.get(this.hammerColIdx));
            String repId = line.get(this.expectColIdx);
            // Insure the test genome is associated with its representative.
            this.expectMap.put(testGenome, repId);
            // Count the hit.
            counters.count(hitGenome);
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
        // This is the number of test genomes with a bad match.
        int badMatchCount = 0;
        // Write the output header.
        writer.println("test_genome_id\ttest_genome_name\ttotal_hits\tmatch_genome_id\tmatch_hits\tmatch_dist\trep_id\trep_hits\trep_dist\tclose\tgood\tbetter");
        // Loop through the test genomes.
        int gCount = 0;
        int gTotal = this.testGenomes.size();
        for (Genome testGenome : this.testGenomes) {
            gCount++;
            log.info("Processing genome {} of {}: {}.", gCount, gTotal, testGenome);
            String testGenomeId = testGenome.getId();
            String firstFields = testGenomeId + "\t" + testGenome.getName();
            // Get the best match.
            var counters = this.hitCountMap.get(testGenomeId);
            var best = counters.getBestEntry();
            if (best == null) {
                // Here there were no hits.  This almost NEVER happens.
                writer.println(firstFields + "0\t\t\t\t\t\t\t\t\t");
                noMatchCount++;
            } else {
                // We will store the rating flags here.
                String closeFlag = "";
                String goodFlag = "";
                String betterFlag = "";
                // Get the genome IDs and hit counts.
                String matchGenomeId = best.getKey();
                int matchHits = best.getCount();
                String repGenomeId = this.expectMap.get(testGenomeId);
                int repHits = counters.getCount(repGenomeId);
                int totalHits = counters.getTotal();
                // Now we compute the distances.  There is always a match distance.
                Measurer testAnalysis = this.distanceComputer.getMeasurer(testGenome);
                Genome matchGenome = this.repGenomes.getGenome(matchGenomeId);
                log.info("Computing distance to best match {}.", matchGenome);
                double matchDist = this.distanceComputer.getDistance(testAnalysis, matchGenome);
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
                    // We need the distance to the expected genome.
                    Genome expectGenome = this.repGenomes.getGenome(repGenomeId);
                    log.info("Computing distance to expected match {}.", expectGenome);
                    repDist = this.distanceComputer.getDistance(testAnalysis, expectGenome);
                    if (repDist >= matchDist) {
                        betterCount++;
                        betterFlag = "Y";
                    } else if (closeFlag.isEmpty()) {
                        // We matched the wrong genome and we are not close, so we count this as a hard miss.
                        badMatchCount++;
                    }
                }
                writer.println(firstFields + String.format("\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                        totalHits, matchGenomeId, matchHits, matchDist, repGenomeId, repHits, repDist,
                        closeFlag, goodFlag, betterFlag));
            }
        }
        // Write the final stats.
        log.info("{} test genomes processed.  {} good hits, {} close hits, {} better hits, {} with no hits, {} hard misses.",
                gTotal, goodCount, closeCount, betterCount, noMatchCount, badMatchCount);
    }

}
