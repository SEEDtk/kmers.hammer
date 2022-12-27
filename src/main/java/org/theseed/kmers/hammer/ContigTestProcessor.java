/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.java.erdb.DbConnection;
import org.theseed.locations.Location;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This command analyzes a single FASTA file to look for hammer hits.  For each hit, it outputs the contig ID and location,
 * the ID of the feature from which the hammer was harvested, and the contig comment.  This information can be used to
 * determine how good the hits were.
 *
 * The positional parameter is the name of the FASTA file.  It is recommended that "genome.bins synth" be used to generate
 * the input file. The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries (default 200)
 * -B	optimal number of kilobases for each sequence group (default 1000)
 *
 * --hType		type of hammer database (default MEMORY)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type (default SQLITE)
 * --headers	comma-delimited list of headers to use for the sequence comment (default "genome_name,rep_id,distance")
 *
 * @author Bruce Parrello
 *
 */
public class ContigTestProcessor extends BaseReportProcessor implements HammerDb.IParms, DbConnection.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ContigTestProcessor.class);
    /** hammer database */
    private HammerDb hammers;
    /** optimal number of base pairs per sequence group */
    private int seqBatchSize;
    /** number of sequences read */
    private int seqsIn;
    /** number of hits output */
    private int hitsOut;
    /** additional headers for sequence comment in output report */
    private String[] headers;

    // COMMAND-LINE OPTIONS

    /** batch size for database or web queries */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for queries")
    private int batchSize;

    /** batch size for sequences */
    @Option(name = "--dnaBatch", aliases = { "-B" }, metaVar = "5", usage = "optimal number of kilobases for sequence groups")
    private int seqBatchKSize;

    /** type of hammer database */
    @Option(name = "--hType", usage = "type of hammer database")
    private HammerDb.Type hammerType;

    /** database engine type */
    @Option(name = "--type", usage = "type of database engine")
    private DbConnection.Type dbEngine;

    /** name of file containing the database */
    @Option(name = "--file", aliases = { "--dbFile", "--dbfile" }, metaVar = "sqlite.db",
            usage = "name of the database file (for SQLITE or flat file)")
    private File dbFile;

    /** database URL for network-based databases */
    @Option(name = "--url", metaVar = "localhost/mainDB", usage = "database resource location (for MySQL)")
    private String dbUrl;

    /** database parameter string */
    @Option(name = "--parms", metaVar = "user=xxx&pass=YYY", usage = "database parameter string (for MySQL")
    private String dbParms;

    /** comma-delimited list of headers for fields in sequence comments (if there are internal tabs) */
    @Option(name = "--headers", metaVar = "genome_id,name", usage = "comma-delimited list of headers to use for sequence comment fields")
    private String headerList;

    /** input FASTA file name */
    @Argument(index = 0, metaVar = "input.fa", usage = "name of the FASTA file containing the sequences to search", required = true)
    private File inFile;

    @Override
    protected void setReporterDefaults() {
        this.batchSize = 200;
        this.hammerType = HammerDb.Type.MEMORY;
        this.dbEngine = DbConnection.Type.SQLITE;
        this.dbFile = null;
        this.dbUrl = null;
        this.dbParms = null;
        this.seqBatchKSize = 1000;
        this.headerList = "genome_name,rep_id,distance";
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify the tuning parameters.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1.");
        if (this.seqBatchKSize < 1)
            throw new ParseFailureException("Dna batch size must be at least 1.");
        // Compute the batch size in base pairs.
        this.seqBatchSize = this.seqBatchKSize * 1024;
        // Verify the database file, if specified.
        if (this.dbFile != null && ! this.dbFile.canRead())
            throw new FileNotFoundException("Database file " + this.dbFile + " is not found or unreadable.");
        // Verify the input file.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        // Parse the comment headers.
        this.headers = StringUtils.split(this.headerList, ',');
        // Set up the hammer database.  We must convert SQL exceptions for compatability.
        try {
            this.hammers = this.hammerType.create(this);
        } catch (SQLException e) {
            throw new IOException("Database error: " + e.toString());
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // This will hold the current sequence batch.  It maps sequence labels to sequences.
        Map<String, Sequence> batch = new HashMap<String, Sequence>(this.seqBatchSize * 4 / 9000 + 10);
        int batchLen = 0;
        // Clear the counters.
        this.seqsIn = 0;
        this.hitsOut = 0;
        int batchCount = 0;
        // Write the output headers.
        writer.println("location\thammer_fid\tworth\t" + StringUtils.join(this.headers, '\t'));
        // Loop through the sequence file.
        try (FastaInputStream inStream = new FastaInputStream(this.inFile)) {
            // Loop through the sequences.
            for (Sequence seq : inStream) {
                // Insure there is room for this sequence.
                if (batchLen >= this.seqBatchSize) {
                    batchCount++;
                    log.info("Processing batch {} of length {}.", batchCount, batchLen);
                    this.processBatch(batch, writer);
                    batch.clear();
                    batchLen = 0;
                }
                // Add the sequence to the batch.
                batch.put(seq.getLabel(), seq);
                batchLen += seq.length();
                this.seqsIn++;
            }
            // Process the residual (if any).
            if (batchLen > 0) {
                log.info("Processing residual batch of length {}.", batchLen);
                this.processBatch(batch, writer);
            }
        }
        log.info("{} sequences read and {} hits found.", this.seqsIn, this.hitsOut);
    }

    /**
     * Process a sequence batch and write out the hits.
     *
     * @param batch		map of sequences to process, keyed by label
     * @param writer	output file for hits
     */
    private void processBatch(Map<String, Sequence> batch, PrintWriter writer) {
        // Get the hits for this batch.
        var hits = this.hammers.findHits(batch.values());
        int hitCount = 0;
        for (HammerDb.Hit hit : hits) {
            final Location hitLoc = hit.getLoc();
            String locString = hitLoc.toSeedString();
            String seqId = hitLoc.getContigId();
            double worth = hit.getStrength();
            String comment = batch.get(seqId).getComment();
            writer.println(locString + "\t" + hit.getFid() + "\t" + worth + "\t" + comment);
            hitCount++;
        }
        writer.flush();
        log.info("{} sequences in batch with {} hits found. {} total sequences read.", batch.size(), hitCount, this.seqsIn);
        this.hitsOut += hitCount;
    }

    @Override
    public File getDbFile() {
        return this.dbFile;
    }

    @Override
    public String getParms() {
        return this.dbParms;
    }

    @Override
    public String getDbUrl() {
        return this.dbUrl;
    }

    @Override
    public DbConnection getConnection()  {
        DbConnection retVal = null;
        try {
            retVal = this.dbEngine.create(this);
        } catch (SQLException e) {
            throw new RuntimeException("Error creating database: " + e.toString());
        }
        return retVal;
    }

    @Override
    public int getBatchSize() {
        return this.batchSize;
    }

}
