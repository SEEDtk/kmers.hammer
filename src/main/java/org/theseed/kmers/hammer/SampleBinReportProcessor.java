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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.fastq.FastqSampleGroup;
import org.theseed.sequence.fastq.ReadStream;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.utils.BaseHammerUsageProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command uses hammers to produce a standard bin report for a set of samples.  The input should be a FASTQ sample group.
 * Currently, this can be a FASTQ sample master directory, a QZA file, or a master directory prepared for binning; however,
 * additional directory structures may be added in the future.
 *
 * The positional parameters are the name of the sample group directory or file and the name of a repgen stats file
 * containing the representative genome names.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 * -t	type of sample group (default FASTQ)
 *
 * --hType		type of hammer database (default MEMORY)
 * --method		voting method to use (default COUNT)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type (default MEMORY)
 * --para		use parallel processing
 * --min		minimum score for an acceptable genome hit (default 10)
 * --qual		minimum acceptable quality for an acceptable sequence (default 10.0)
 * --seqBatch	maximum number of kilobases to process at one time (default 1000)
 *
 * @author Bruce Parrello
 *
 */
public class SampleBinReportProcessor extends BaseHammerUsageProcessor {

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
    /** number of sequences read */
    private int seqsIn;
    /** number of sequences rejected due to quality */
    private int qualReject;
    /** number of batches processed */
    private int batchCount;
    /** number of genomes rejected due to low score */
    private int badScores;
    /** number of genomes output for all samples */
    private int goodScores;
    /** maximum batch size in base pairs */
    private int maxBatchDnaSize;

    // COMMAND-LINE OPTIONS

    /** type of input sample group */
    @Option(name = "--sType", aliases = { "-t" }, usage = "type of input sample group")
    private FastqSampleGroup.Type groupType;

    /** if specified, parallel processing will be used */
    @Option(name = "--para", usage = "if specified, samples will be processed in parallel")
    private boolean paraFlag;

    /** minimum acceptable hammer score */
    @Option(name = "--min", metaVar = "20.0", usage = "minimum acceptable score for a genome presence")
    private double minScore;

    /** minimum acceptable quality for an input sequence */
    @Option(name = "--qual", metaVar = "30.0", usage = "minimum acceptable sequence quality (0 to 99)")
    private double minQual;

    /** maximum number of kilobases */
    @Option(name = "--seqBatch", metaVar = "2000", usage = "maximum number of kilobases to process in a batch")
    private int seqBatchSize;

    /** input sample group directory (or file) */
    @Argument(index = 0, metaVar = "inDir", usage = "input sample group directory (or file for certain types)",
            required = true)
    private File inDir;

    /** repgen stats file containing genome names */
    @Argument(index = 1, metaVar = "repgen.stats.tbl", usage = "repgen stats file containing representative genome names")
    private File repStatsFile;

    @Override
    protected void setHammerDefaults() {
        this.paraFlag = false;
        this.groupType = FastqSampleGroup.Type.FASTQ;
        this.minScore = 10.0;
        this.minQual = 10.0;
        this.seqBatchSize = 1000;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Validate the minimum score.
        if (this.minScore < 0.0)
            throw new ParseFailureException("Minimum score cannot be negative.");
        // Insure the sample group directory exists.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input sample group " + this.inDir + " not found in file system.");
        // Insure it is valid.
        FileFilter filter = this.groupType.getFilter();
        if (! filter.accept(this.inDir))
            throw new IOException(this.inDir + " does not appear to be a valid sample group of type " + this.groupType + ".");
        // Validate the sequence batch size.
        if (this.seqBatchSize < 1)
            throw new ParseFailureException("Maximum sequence batch size must be at least 1.");
        this.maxBatchDnaSize = this.seqBatchSize * 1024;
        // Read in the genome names.
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + " is not found or unreadable.");
        this.genomeMap = TabbedLineReader.readMap(this.repStatsFile, "1", "2");
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
        this.goodScores = 0;
        // Write the report header.
        writer.println("sample_id\trepgen_id\trepgen_name\tcount");
        // Start by connecting to the sample group.
        try (FastqSampleGroup sampleGroup = this.groupType.create(this.inDir)) {
            Set<String> samples = sampleGroup.getSamples();
            log.info("{} samples found in {}.", samples.size(), this.inDir);
            // Now we loop through the samples.
            Stream<String> sampleStream = samples.stream();
            if (this.paraFlag)
                sampleStream = sampleStream.parallel();
            log.info("Processing samples.");
            sampleStream.forEach(x -> this.processSample(sampleGroup, x, writer));
        }
        if (log.isInfoEnabled() && this.sampleCount > 0) {
            Duration allTime = Duration.ofMillis(this.processTime);
            Duration perSample = Duration.ofMillis(this.processTime / this.sampleCount);
            log.info("Process time was {}.  {} per sample.", allTime, perSample);
            log.info("{} sequences read, {} rejected, {} batches processed.", this.seqsIn, this.qualReject, this.batchCount);
            log.info("{} sample/genome pairs output, {} rejected due to low score.", this.goodScores, this.badScores);
        }
    }

    /**
     * Scan a sample for hammers and output the genome scores.
     *
     * @param sampleId	ID of the sample to scan
     * @param writer	output writer for the results
     */
    private void processSample(FastqSampleGroup sampleGroup, String sampleId, PrintWriter writer) {
        long start = System.currentTimeMillis();
        // Get the read stream.
        ReadStream inStream = this.getStream(sampleGroup, sampleId);
        // This timer is used to space out log messages.
        long lastMessage = 0;
        // The hammer database operates on a collection of sequences, while we have a stream of
        // sequence reads.  In general, the reads will be very long for FASTA samples, with varying
        // coverage, and very short for FASTQ samples, with uniform coverage.  To handle both
        // these cases efficiently, we collect the reads into batches with identical coverage.  For
        // FASTA, these batches will generally be size 1, but with lots of DNA.  For FASTQ, the
        // batches will be large size, so also with lots of DNA.
        int batchSize = 0;
        List<Sequence> batch = new ArrayList<Sequence>(batchSize);
        double batchCoverage = 1.0;
        // Start with an empty weight map.
        WeightMap results = new WeightMap(this.genomeMap.size() * 4 / 3 + 1);
        // Get local versions of the counters.
        int mySeqsIn = 0;
        int myBatchCount = 0;
        int myQualReject = 0;
        // Loop through the stream.
        while (inStream.hasNext()) {
            SeqRead seqRead = inStream.next();
            mySeqsIn++;
            // Filter on quality.
            if (seqRead.getQual() < this.minQual)
                myQualReject++;
            else {
                double coverage = seqRead.getCoverage();
                // Insure there is room in this batch for this sequence.
                if (batchSize >= this.maxBatchDnaSize || coverage != batchCoverage) {
                    // Here we must process the current batch.  Get the scores and
                    // add them in.
                    myBatchCount++;
                    WeightMap scores = this.processBatch(batch);
                    results.accumulate(scores, batchCoverage);
                    // Set up for the new sequence.
                    batchCoverage = coverage;
                    batch.clear();
                    if (log.isInfoEnabled() && System.currentTimeMillis() - lastMessage > 5000) {
                        lastMessage = System.currentTimeMillis();
                        double rate = mySeqsIn * 1000.0 / (lastMessage - start);
                        log.info("{} sequences read, {} rejected from {}, {} sequences/second.", mySeqsIn, myQualReject, sampleId, rate);
                    }
                }
                // Convert the read to a sequence and add it to the batch.
                Sequence seq = new Sequence(seqRead.getLabel(), "", seqRead.getSequence());
                batch.add(seq);
                batchSize += seq.length();
            }
        }
        // Process the residual batch.
        if (batch.size() > 0) {
            myBatchCount++;
            WeightMap scores = this.processBatch(batch);
            results.accumulate(scores, batchCoverage);
        }
        log.info("Sample {} contained {} sequences in {} batches.  {} were rejected due to quality.  {} genomes scored.",
                sampleId, mySeqsIn, myBatchCount, myQualReject, results.size());
        // Output the genome weights.
        this.writeSample(writer, sampleId, results);
        // Update the counts.
        synchronized(this) {
            this.batchCount += myBatchCount;
            this.qualReject += myQualReject;
            this.seqsIn += mySeqsIn;
            this.processTime += System.currentTimeMillis() - start;
            this.sampleCount++;
        }
    }

    /**
     * @return a read strean for the specified sample
     *
     * @param sampleGroup	master sample group
     * @param sampleId		ID of the desired sample
     */
    private ReadStream getStream(FastqSampleGroup sampleGroup, String sampleId) {
        ReadStream retVal;
        try {
            retVal = sampleGroup.sampleIter(sampleId);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        return retVal;
    }

    /**
     * Process a batch of sequences and return the genome weights.
     *
     * @param batch		 batch of sequences to process
     *
     * @return a weight map containing the score for each genome in this batch.
     */
    private WeightMap processBatch(List<Sequence> batch) {
        return this.hammers.findClosest(batch);
    }

    /**
     * Write the specified sample to the output report.
     *
     * @param writer	output writer for the report
     * @param sampleId	ID of the sample being processed
     * @param results	scores for this sample
     */
    private synchronized void writeSample(PrintWriter writer, String sampleId, WeightMap results) {
        log.info("Writing results for sample {}.", sampleId);
        var scores = results.sortedCounts();
        for (var score : scores) {
            double scoreVal = score.getCount();
            if (scoreVal < this.minScore)
                this.badScores++;
            else {
                this.goodScores++;
                String genomeId = score.getKey();
                String genomeName = this.genomeMap.getOrDefault(genomeId, "<unknown>");
                writer.println(sampleId + "\t" + genomeId + "\t" + genomeName + "\t" + score.getCount());
            }
        }
        // Try to keep the output aligned on sample boundaries so we can restart.
        writer.flush();
    }


}