/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.ForkJoinPool;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.ClassStrategy;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.fastq.FastqSampleGroup;
import org.theseed.sequence.fastq.ReadStream;
import org.theseed.sequence.fastq.SampleDescriptor;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.stats.WeightMap;
import org.theseed.proteins.hammer.ScoreMap;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command uses hammers to produce a standard bin report for a set of
 * samples. The input should be a FASTQ sample group. Currently, this can be a
 * FASTQ sample master directory, a QZA file, or a master directory prepared for
 * binning; however, additional directory structures may be added in the future.
 *
 * The positional parameters are the name of the sample group directory or file
 * and the name of a repgen stats file containing the representative genome
 * names. The command-line options are as follows.
 *
 *  -h	display command-line usage
 *  -v 	display more frequent log messages -o output file (if not STDOUT)
 *  -b	batch size for queries
 *  -t	type of sample group (default FASTQ)
 *
 *  --hType		type of hammer database (default MEMORY)
 *  --method	voting method to use (default COUNT)
 *  --strategy	scoring strategy to use (default HITS)
 *  --file		file containing hammer database (either SQLite database or hammer flat file)
 *  --url		URL of database (host and name, MySQL only)
 *  --parms 	database connection parameter string (MySQL only)
 *  --type		database engine type (default MEMORY)
 *  --para 		number of parallel threads to use
 *  --min 		minimum score for an acceptable role/genome hit (default 10)
 *  --qual 		minimum quality for an acceptable sequence hit (default 0.7)
 *  --seqBatch 	maximum number of kilobases to process at one time (default 1000)
 *  --diff 		minimum number of additional hits for a region-based classification
 *  --resume 	name of an existing report file containing data from samples already run
 *  --unscaled 	if specified, hit weights will not be scaled by length and coverage
 *  --filter	maximum number of expected bad base pairs for an acceptable read
 *  --minG		minimum score for an acceptable genome hit (default 1.0)
 *
 * @author Bruce Parrello
 *
 */
public class SampleBinReportProcessor extends BaseHammerUsageProcessor implements ClassStrategy.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SampleBinReportProcessor.class);
    /** hammer database */
    private HammerDb hammers;
    /** map of genome IDs to names */
    private Map<String, String> genomeMap;
    /** accumulator for time spent scanning for hammers */
    private long processTime;
    /** number of samples processed */
    private int sampleCount;
    /** set of bad samples */
    private Set<String> badSamples;
    /** number of sequences read */
    private long seqsIn;
    /** number of sequences rejected due to quality */
    private long qualReject;
    /** number of batches processed */
    private long batchCount;
    /** number of role/genome pairs rejected due to low score */
    private long badScores;
    /** number of role/genome pairs output for all samples */
    private long goodScores;
    /** number of genome presences output */
    private long goodGenomes;
    /** number of genome presences rejected due to low score */
    private long badGenomes;
    /** maximum batch size in base pairs */
    private int maxBatchDnaSize;
    /** classification strategy helper */
    private ClassStrategy strategy;
    /** weight map size to use */
    private int mapSize;
    /** set of roles being tracked, in order */
    private List<String> roleSet;
    /** custom thread pool for parallel processing */
    private ForkJoinPool threadPool;

    // COMMAND-LINE OPTIONS

    /** type of input sample group */
    @Option(name = "--sType", aliases = { "-t" }, usage = "type of input sample group")
    private FastqSampleGroup.Type groupType;

    /** type of classification strategy */
    @Option(name = "--strategy", usage = "strategy for classifying sequences")
    private ClassStrategy.Type strategyType;

    /** number of threads to use in parallel processing */
    @Option(name = "--para", metaVar = "60", usage = "maximum number of threads to run in parallel")
    private int maxThreads;

    /** minimum acceptable hammer score for role in genome */
    @Option(name = "--min", metaVar = "20.0", usage = "minimum acceptable score for a role/genome presence")
    private double minScore;

    /** minimum acceptable hammer score for whole genome */
    @Option(name = "--minG", metaVar = "10.0", usage = "minimum acceptable score for a genome presence")
    private double minGenomeScore;

    /** minimum number of additional hits for a region-based determination to count */
    @Option(name = "--diff", metaVar = "5", usage = "minimum number of additional hits required to qualify as a regional group classification")
    private double minDiff;

    /** maximum number of expected errors allowed in an acceptable read */
    @Option(name = "--filter", metaVar = "2.0", usage = "number of expected errors at or beyond which a read is rejected")
    private double badBaseFilter;

    /** minimum acceptable quality for an sequence hit */
    @Option(name = "--qual", metaVar = "0.9", usage = "minimum acceptable sequence hit quality (0 to 1)")
    private double minQual;

    /** maximum number of kilobases per batch */
    @Option(name = "--seqBatch", metaVar = "2000", usage = "maximum number of kilobases to process in a batch")
    private int seqBatchSize;

    /** name of previous report file */
    @Option(name = "--resume", metaVar = "binReport.old.tbl", usage = "if specified, the name of a previous report containing results for samples that do not need to be rerun (cannot also be the output file)")
    private File resumeFile;

    /** scale-suppression flag */
    @Option(name = "--unscaled", usage = "if specified, hit weights will NOT be scaled by length and coverage")
    private boolean unScaled;

    /** input sample group directory (or file) */
    @Argument(index = 0, metaVar = "inDir", usage = "input sample group directory (or file for certain types)", required = true)
    private File inDir;

    /** repgen stats file containing genome names */
    @Argument(index = 1, metaVar = "repgen.stats.tbl", usage = "repgen stats file containing representative genome names", required = true)
    private File repStatsFile;

    /**
     * This is a small utility object that tracks hit/read statistics for a sample.
     */
    protected class SampleStats {

        /** number of sequences with hits */
        private int hitSeqCount;
        /** number of hits */
        private int hitCount;
        /** number of ambiguous reads */
        private int ambigCount;

        protected SampleStats() {
            this.hitSeqCount = 0;
            this.hitCount = 0;
            this.ambigCount = 0;
        }

    }

    @Override
    protected void setHammerDefaults() {
        this.groupType = FastqSampleGroup.Type.FASTQ;
        this.minScore = 10.0;
        this.minQual = 0.7;
        this.seqBatchSize = 1000;
        this.strategyType = ClassStrategy.Type.HITS;
        this.minDiff = 2;
        this.minGenomeScore = 1.0;
        this.resumeFile = null;
        this.badBaseFilter = 1.0;
        this.maxThreads = Runtime.getRuntime().availableProcessors();
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Validate the minimum scores.
        if (this.minScore < 0.0)
            throw new ParseFailureException("Minimum role-score cannot be negative.");
        if (this.minGenomeScore < 0.0)
            throw new ParseFailureException("Minimum genome-score cannot be negative.");
        // Validate the bad-base filter.
        if (this.badBaseFilter <= 0.0)
            throw new ParseFailureException("Bad-base filter limit must be positive.");
        // Validate the minimum difference.
        if (this.minDiff <= 0.0)
            throw new ParseFailureException("Minimum hit difference must be greater than 0.");
        // Validate the minimum quality for a hit.
        if (this.minQual <= 0.0 || this.minQual >= 1.0)
            throw new ParseFailureException("Minimum quality must be between 0 and 1.");
        // Insure the sample group directory exists.
        if (!this.inDir.exists())
            throw new FileNotFoundException("Input sample group " + this.inDir + " not found in file system.");
        // Insure it is valid.
        FileFilter filter = this.groupType.getFilter();
        if (!filter.accept(this.inDir))
            throw new IOException(
                    this.inDir + " does not appear to be a valid sample group of type " + this.groupType + ".");
        // Validate the sequence batch size.
        if (this.seqBatchSize < 1)
            throw new ParseFailureException("Maximum sequence batch size must be at least 1.");
        this.maxBatchDnaSize = this.seqBatchSize * 1024;
        // Read in the genome names.
        if (!this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + " is not found or unreadable.");
        this.genomeMap = TabbedLineReader.readMap(this.repStatsFile, "1", "2");
        // Check for a resume file.
        if (this.resumeFile != null && !this.resumeFile.canRead())
            throw new FileNotFoundException("Resume file " + this.resumeFile + " is not found or unreadable.");
        // Compute the weight map size to use.
        this.mapSize = this.genomeMap.size() * 4 / 3 + 1;
        // Create the strategy helper.
        this.strategy = this.strategyType.create(this);
        int maxCores = Runtime.getRuntime().availableProcessors();
        if (this.maxThreads > 1 && this.maxThreads > maxCores) {
            log.warn("Too many threads specified:  reducing from {} to {}.", this.maxThreads, maxCores);
            this.maxThreads = maxCores;
        }
        // Create the custom thread pool.
        if (this.maxThreads == 1)
            this.threadPool = null;
        else {
            this.threadPool = new ForkJoinPool(this.maxThreads);
            log.info("Parallel processing selected with {} threads.", this.maxThreads);
        }
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        this.hammers = hammerDb;
        // Clear the tracking fields.
        this.processTime = 0;
        this.sampleCount = 0;
        this.seqsIn = 0;
        this.batchCount = 0;
        this.qualReject = 0;
        this.badScores = 0;
        this.badGenomes = 0;
        this.goodGenomes = 0;
        this.goodScores = 0;
        this.badSamples = new TreeSet<String>();
        // Get the set of roles.
        this.roleSet = new ArrayList<String>(this.hammers.getRoles());
        // Write the report header.
        writer.println("sample_id\trepgen_id\trepgen_name\tcount\troles\t" + StringUtils.join(this.roleSet, '\t'));
        // Start by connecting to the sample group.
        try (FastqSampleGroup sampleGroup = this.groupType.create(this.inDir)) {
            Set<String> samples = sampleGroup.getSampleIDs();
            log.info("{} samples found in {}.", samples.size(), this.inDir);
            // Now process the resume file (if any).
            if (this.resumeFile != null)
                this.processResume(samples, writer);
            // Now we loop through the samples.
            log.info("Processing samples.");
            // Are we parallelizing?
            if (this.threadPool == null) {
                // No, use a normal stream.
                sampleGroup.stream(samples, false).forEach(x -> this.processSample(x, writer));
            } else try {
                // Yes, use a parallel stream.
                this.threadPool.submit(() -> sampleGroup.stream(samples, true).forEach(x -> this.processSample(x, writer))).get();
            } finally {
                this.threadPool.shutdown();
            }
        }
        if (log.isInfoEnabled() && this.sampleCount > 0) {
            Duration allTime = Duration.ofMillis(this.processTime);
            Duration perSample = Duration.ofMillis(this.processTime / this.sampleCount);
            log.info("Process time was {}.  {} per sample.", allTime, perSample);
            log.info("{} sequences read, {} rejected, {} batches processed.", this.seqsIn, this.qualReject,
                    this.batchCount);
            log.info("{} role/genome pairs output, {} role hits rejected due to low score, "
                    + "{} genome hits output, {} rejected due to low score, {} bad-quality hits ignored.",
                    this.goodScores, this.badScores, this.goodGenomes, this.badGenomes,
                    this.hammers.getBadQual());
        }
        if (!this.badSamples.isEmpty())
            log.warn("The following samples had no data: {}", StringUtils.join(badSamples, ", "));
    }

    /**
     * Process the resume file. We copy the file to the output and remove the sample
     * IDs found from the list of samples to process.
     *
     * @param samples list of samples to process (will be modified)
     * @param writer  print writer for the output report
     *
     * @throws IOException
     */
    private void processResume(Set<String> samples, PrintWriter writer) throws IOException {
        try (TabbedLineReader resumeStream = new TabbedLineReader(this.resumeFile)) {
            if (resumeStream.findColumn("sample_id") != 0)
                throw new IOException(this.resumeFile + " does not appear to be a bin report file.");
            // Loop through the resume file, copying report data.
            log.info("Restoring old sample data from {}.", this.resumeFile);
            int lineIn = 0;
            Set<String> newSamples = new HashSet<String>(samples.size() * 4 / 3 + 1);
            for (var line : resumeStream) {
                String newSample = line.get(0);
                newSamples.add(newSample);
                writer.println(line.toString());
                lineIn++;
            }
            log.info("{} lines read containing {} samples.", lineIn, newSamples.size());
            // Remove the samples found here from the set of samples to process.
            samples.removeAll(newSamples);
            log.info("{} samples remaining after resume processing.", samples.size());
        }
    }

    /**
     * Scan a sample for hammers and output the genome scores.
     *
     * @param sample descriptor for the sample to scan
     * @param writer output writer for the results
     */
    private void processSample(SampleDescriptor sample, PrintWriter writer) {
        long start = System.currentTimeMillis();
        // Get the sample ID.
        String sampleId = sample.getId();
        log.info("Processing sample {}.", sampleId);
        // The hammer database operates on a collection of sequences, while we have a stream of
        // sequence reads. In general, the reads will be very long for FASTA samples, with varying
        // coverage, and very short for FASTQ samples, with uniform coverage. To handle both
        // these cases efficiently, we collect the reads into batches with identical coverage. For
        // FASTA, these batches will generally be size 1, but with lots of DNA. For FASTQ, the
        // batches will be large size, so also with lots of DNA.
        int batchSize = 0;
        var batch = new HashMap<String, SeqRead.Part>();
        double batchCoverage = 1.0;
        // Start with an empty weight map.
        ScoreMap results = new ScoreMap(this.mapSize);
        // Get local versions of the counters.
        int mySeqsIn = 0;
        int myBatchCount = 0;
        int rejectCount = 0;
        // Denote the sample is bad.
        boolean badSample = true;
        // Create a stats object to pass to the batch processor.
        SampleStats stats = new SampleStats();
        // Get the read stream.
        try (ReadStream inStream = this.getStream(sample)) {
            // This timer is used to space out log messages.
            long lastMessage = 0;
            // Loop through the stream.
            while (inStream.hasNext()) {
                SeqRead seqRead = inStream.next();
                mySeqsIn++;
                double coverage = seqRead.getCoverage();
                double expError = seqRead.getExpectedErrors();
                if (expError >= this.badBaseFilter) {
                    rejectCount++;
                    log.debug("{} rejected due to expected error {}.", seqRead.getLabel(), expError);
                } else {
                    // Insure there is room in this batch for this sequence.
                    if (batchSize >= this.maxBatchDnaSize || coverage != batchCoverage) {
                        // Here we must process the current batch. Get the scores and
                        // add them in.
                        myBatchCount++;
                        ScoreMap scores = this.processBatch(batch, batchCoverage, stats);
                        results.accumulate(scores);
                        // Set up for the new sequence.
                        batchCoverage = coverage;
                        batch.clear();
                        batchSize = 0;
                        if (log.isInfoEnabled() && System.currentTimeMillis() - lastMessage > 10000) {
                            lastMessage = System.currentTimeMillis();
                            double rate = mySeqsIn * 1000.0 / (lastMessage - start);
                            log.info("{} sequences read in {}, {} sequences/second.", mySeqsIn, sampleId, rate);
                        }
                    }
                    // Convert the read to a sequence and add it to the batch.
                    final String seqId = seqRead.getLabel();
                    SeqRead.Part seq = seqRead.getSequence();
                    batch.put(seqId, seq);
                    batchSize += seq.length();
                }
            }
            // Process the residual batch.
            if (batch.size() > 0) {
                myBatchCount++;
                ScoreMap scores = this.processBatch(batch, batchCoverage, stats);
                results.accumulate(scores);
            }
            // No errors, so the sample is good.
            badSample = false;
        } catch (UncheckedIOException e) {
            IOException cause = e.getCause();
            String message;
            if (cause != null)
                message = cause.toString();
            else
                message = e.getMessage();
            log.error("ERROR processing sample {}: {}", sampleId, message);
        }
        if (badSample)
            log.error("Sample {} suppressed due to error.", sampleId);
        else if (mySeqsIn == 0) {
            log.warn("WARNING: sample {} had no data.", sampleId);
            badSample = true;
        } else if (stats.hitSeqCount == 0) {
            log.warn("WARNING: sample {} had no classifiable sequences.", sampleId);
            badSample = true;
        } else {
            log.info(
                    "Sample complete: {} contained {} sequences in {} batches.  {} genomes scored.  {} ambiguous, {} rejected, {} sequences classified, {} per sequence.",
                    sampleId, mySeqsIn, myBatchCount, results.size(), stats.ambigCount, rejectCount, stats.hitSeqCount,
                    ((double) stats.hitCount) / stats.hitSeqCount);
            // Output the genome weights.
            this.writeSample(writer, sampleId, results);
        }
        // Update the counts.
        synchronized (this) {
            this.batchCount += myBatchCount;
            this.seqsIn += mySeqsIn;
            this.processTime += System.currentTimeMillis() - start;
            this.sampleCount++;
            if (badSample)
                this.badSamples.add(sampleId);
            this.qualReject += rejectCount;
        }
    }

    /**
     * @return a read strean for the specified sample
     *
     * @param sample sample to read
     */
    private ReadStream getStream(SampleDescriptor sample) {
        ReadStream retVal;
        try {
            retVal = sample.reader();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        return retVal;
    }

    /**
     * Process a batch of sequences and return the genome weights.
     *
     * @param batch    batch of sequences to process
     * @param coverage coverage for all the sequences in the batch
     *
     * @return a weight map containing the score for each genome in this batch's
     *         sequences
     */
    private ScoreMap processBatch(Map<String, SeqRead.Part> batch, double coverage, SampleStats stats) {
        ScoreMap retVal = new ScoreMap(this.mapSize);
        // Get all the hits for the batch, sorted by location.
        SortedSet<HammerDb.Hit> hits = this.hammers.findHits(batch.values(), this.minQual);
        // Loop through the hits. For each sequence, we compute its scores and
        // accumulate them in the return map. This requires batching the hits. Because of the sort, all
        // the hits for a sequence will be together. We need to remember the ID and length of the
        // current sequence.
        String seqId = "";
        int seqLen = 0;
        Collection<HammerDb.Hit> hitSet = new ArrayList<HammerDb.Hit>(hits.size());
        for (var hit : hits) {
            String hitSeqId = hit.getLoc().getContigId();
            if (! hitSeqId.contentEquals(seqId)) {
                // This hit is for a new batch. Process the old one.
                if (hitSet.size() > 0) {
                    this.updateMap(coverage, stats, retVal, seqLen, hitSet);
                    hitSet.clear();
                }
                seqId = hitSeqId;
                seqLen = batch.get(hitSeqId).length();
            }
            // Add this hit to the current batch.
            hitSet.add(hit);
        }
        // Check for a residual batch.
        if (hitSet.size() > 0) {
            this.updateMap(coverage, stats, retVal, seqLen, hitSet);
        }
        // Return the accumulated scores.
        return retVal;
    }

    /**
     * Update the master weight map from the hits for a sequence.
     *
     * @param coverage  sequence coverage
     * @param stats     statistics for this sample
     * @param masterMap master weight map to update
     * @param seqLen    length of the sequence
     * @param hitSet    set of hits against the sequence
     */
    public void updateMap(double coverage, SampleStats stats, ScoreMap masterMap, int seqLen,
            Collection<HammerDb.Hit> hitSet) {
        ScoreMap batchMap = this.strategy.computeScores(hitSet, seqLen, coverage);
        if (batchMap.size() <= 0)
            stats.ambigCount++;
        else {
            masterMap.accumulate(batchMap);
            stats.hitSeqCount++;
            stats.hitCount += hitSet.size();
        }
    }

    /**
     * Write the specified sample to the output report.
     *
     * @param writer   output writer for the report
     * @param sampleId ID of the sample being processed
     * @param results  scores for this sample
     */
    private synchronized void writeSample(PrintWriter writer, String sampleId, ScoreMap results) {
        log.info("Writing results for sample {}.", sampleId);
        var scores = results.sortedCounts();
        for (var score : scores) {
            WeightMap roleScores = score.getRoleScores();
            // We need the total score, the number of roles found, and the output string for the role counts.
            double totalScore = 0.0;
            int roleCount = 0;
            StringBuffer roleString = new StringBuffer(100);
            // Only keep the roles with high enough scores.
            for (String role : this.roleSet) {
                double scoreVal = roleScores.getCount(role);
                if (scoreVal < this.minScore) {
                    this.badScores++;
                    roleString.append("\t0.0");
                } else {
                    this.goodScores++;
                    roleCount++;
                    roleString.append('\t').append(scoreVal);
                    totalScore += scoreVal;
                }
            }
            if (totalScore < this.minGenomeScore) {
                // Here the total hit is too small.
                this.badGenomes++;
            } else {
                // Here we have an output line.  Get the genome data.
                this.goodGenomes++;
                String genomeId = score.getKey();
                String genomeName = this.genomeMap.getOrDefault(genomeId, "<unknown>");
                // Output the hit data.
                writer.println(sampleId + "\t" + genomeId + "\t" + genomeName + "\t" + totalScore + "\t" + roleCount
                        + roleString.toString());
            }
        }
        // Try to keep the output aligned on sample boundaries so we can restart.
        writer.flush();
        log.info("Sample {} completed.", sampleId);
    }

    @Override
    public double getMinDiff() {
        return this.minDiff;
    }

    @Override
    public boolean isScaled() {
        return !this.unScaled;
    }

}
