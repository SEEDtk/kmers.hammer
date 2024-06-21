/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.locations.Location;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.fastq.FastqSampleGroup;
import org.theseed.sequence.fastq.ReadStream;
import org.theseed.sequence.fastq.SampleDescriptor;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command is used to analyze the quality of hammer hits against reads.  The input should be FASTQ samples in a directory,
 * and the output will be a detail report to the standard output listing the individual hits.  A blacklist of genome IDs can be
 * specified for each sample, and hits to these genomes will be skipped on the report.
 *
 * The blacklist is specified in a tab-delimited file with labels.  The first column contains sample IDs, and the second contains a
 * comma-delimited list of genome IDs to be skipped in the detail report.  In essence, the blacklist contains the genome IDs we
 * expect to find in the sample, and other genome IDs are considered suspicious.
 *
 * The positional parameter is the name of the input sample directory.
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
 * --blacklist	blacklist file, containing genome IDs to skip in the detail report for each sample
 * --qual 		minimum quality for an acceptable sequence hit (default 0.7)
 * --filter		maximum number of expected bad base pairs for an acceptable read (default 1.0)
 * --para		run the samples in parallel
 * --good		if specified, an output file to contain the blacklisted hits
 *
 * @author Bruce Parrello
 *
 */
public class ReadTestProcessor extends BaseHammerUsageProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ReadTestProcessor.class);
    /** map of sample IDs to black list sets */
    private Map<String, Set<String>> blackListMap;
    /** hash size for sequence batches */
    private int hashSize;
    /** output stream for blocklisted hits */
    private PrintWriter goodStream;
    /** header line for output files */
    private static final String HEADER_LINE = "sample_id\thammer\thammer_fid\tgenome_id\tread_id\tweight\tquality\trole";

    // COMMAND-LINE OPTIONS

    /** blacklist file name */
    @Option(name = "--blacklist", metaVar = "sampleBlackList.tbl", usage = "optional file containing sample IDs and genomes considered good for each")
    private File blackFile;

    /** TRUE if we should run the samples in parallel */
    @Option(name = "--para", usage = "if specified, parallel-processing will be used to improve performance")
    private boolean paraFlag;

    /** maximum number of expected errors allowed in an acceptable read */
    @Option(name = "--filter", metaVar = "2.0", usage = "number of expected errors at or beyond which a read is rejected")
    private double badBaseFilter;

    /** minimum acceptable quality for a sequence hit */
    @Option(name = "--qual", metaVar = "0.9", usage = "minimum acceptable sequence hit quality (0 to 1)")
    private double minQual;

    /** if specified, the name of an output file for the blacklisted hits */
    @Option(name = "--good", metaVar = "expectedHits.tbl", usage = "optional output file for blacklisted (expected) hits")
    private File goodFile;

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing FASTQ samples", required = true)
    private File inDir;


    @Override
    protected void setHammerDefaults() {
        this.blackFile = null;
        this.paraFlag = false;
        this.badBaseFilter = 1.0;
        this.minQual = 0.7;
        this.goodFile = null;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Insure the input directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input FASTQ directory sample group " + this.inDir + " is not found or is invalid.");
        // Process the black list.
        if (this.blackFile == null) {
            log.info("No blacklist file specified.");
            this.blackListMap = Collections.emptyMap();
        } else if (! this.blackFile.canRead())
            throw new FileNotFoundException("Blacklist file " + this.blackFile + " is not found or unreadable.");
        else {
            // Read in the blacklist.
            Map<String, Set<String>> blackMap = this.readBlackList(this.blackFile);
            this.blackListMap = blackMap;
        }
        // Set up the expected-hits output file.
        if (this.goodFile == null)
            this.goodStream = null;
        else {
            this.goodStream = new PrintWriter(this.goodFile);
            log.info("Blacklisted hits will be written to {}.", this.goodFile);
        }
        // Compute the hash size.
        this.hashSize = this.getBatchSize() * 4 / 3 + 1;
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        try {
            // Write the output header.
            writer.println(HEADER_LINE);
            if (this.goodStream != null)
                this.goodStream.println(HEADER_LINE);
            // Create the sample group.
            try (FastqSampleGroup samples = FastqSampleGroup.Type.FASTQ.create(this.inDir)) {
                log.info("{} samples found in {}.", samples.size());
                samples.stream(this.paraFlag).forEach(x -> this.processSample(hammerDb, x, writer));
            }
        } finally {
            // Insure the good-hits file is closed.
            if (this.goodStream != null)
                this.goodStream.close();
        }
    }

    /**
     * Find all the suspicious hits in a particular sample's reads.
     *
     * @param hammers	hammer database
     * @param sample	sample to examine
     * @param writer	output writer for report
     */
    private void processSample(HammerDb hammers, SampleDescriptor sample, PrintWriter writer) {
        String sampleId = sample.getId();
        log.info("Processing sample {}.", sampleId);
        // Loop through the sample's read stream.
        try (ReadStream readStream = sample.reader()) {
            // This will hold the current batch of reads.
            Map<String, SeqRead.Part> batch = new HashMap<String, SeqRead.Part>(this.hashSize);
            int readCount = 0;
            int batchCount = 0;
            int rejectCount = 0;
            for (SeqRead read : readStream) {
                readCount++;
                if (batch.size() >= this.getBatchSize()) {
                    // Here we have to process a batch.
                    batchCount++;
                    log.info("Processing batch {} for {}, {} total reads.", batchCount, sampleId, readCount);
                    this.processBatch(hammers, batch, sampleId, writer);
                    batch.clear();
                }
                // Check the read for errors.
                if (read.getExpectedErrors() >= this.badBaseFilter)
                    rejectCount++;
                else
                    batch.put(read.getLabel(), read.getSequence());
            }
            // Process the residual.
            if (batch.size() > 0) {
                log.info("Processing residual batch of size {} for {}.", batch.size(), sampleId);
                this.processBatch(hammers, batch, sampleId, writer);
                batchCount++;
            }
            log.info("{} reads and {} batches in sample {}. {} rejected.", readCount, batchCount, sampleId, rejectCount);
        } catch (IOException e) {
            // Convert to untracked so we can use in streaming.
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Process a single batch.  We find all the hits and output the ones we want.
     *
     * @param hammers	hammer database
     * @param batch		batch of sequences to process
     * @param sampleId	ID of the current sample
     * @param writer	output writer for report
     */
    private void processBatch(HammerDb hammers, Map<String, SeqRead.Part> batch, String sampleId, PrintWriter writer) {
        // Get the blacklist.
        Set<String> blacklist = this.blackListMap.getOrDefault(sampleId, Collections.emptySet());
        // Find all the hits for this batch.
        Collection<HammerDb.Hit> hits = hammers.findHits(batch.values(), this.minQual);
        for (HammerDb.Hit hit : hits) {
            // Check the genome ID to insure we want this hit.
            String genomeId = hit.getGenomeId();
            if (! blacklist.contains(genomeId))
                this.writeHit(hammers, batch, sampleId, writer, hit, genomeId);
            else if (this.goodStream != null)
                this.writeHit(hammers, batch, sampleId, this.goodStream, hit, genomeId);
        }
    }

    /**
     * Write the information for a hammer hit to the specified output stream.
     *
     * @param hammers	main hammer database
     * @param batch		current read batch
     * @param sampleId	ID of the sample containing the hit
     * @param writer	output stream for recording the hit
     * @param hit		hit descriptor
     * @param genomeId	repgen ID for the hammer
     */
    public void writeHit(HammerDb hammers, Map<String, SeqRead.Part> batch, String sampleId, PrintWriter writer,
            HammerDb.Hit hit, String genomeId) {
        // Get the basic output fields we need.
        String hammerFid = hit.getFid();
        Location loc = hit.getLoc();
        String readId = loc.getContigId();
        double weight = hit.getStrength();
        String hammer = hit.getHammer();
        // Get the sequence hit.
        SeqRead.Part hitPart = batch.get(readId);
        double qual = SeqRead.qualChance(hitPart.getQual(), loc.getLeft() - 1, hammers.getKmerSize());
        // Write the hit.
        synchronized (writer) {
            writer.println(sampleId + "\t" + hammer + "\t" + hammerFid + "\t" + genomeId + "\t" + readId + "\t"
                    + Double.toString(weight) + "\t" + Double.toString(qual) + "\t" + hit.getRole());
        }
    }

}
