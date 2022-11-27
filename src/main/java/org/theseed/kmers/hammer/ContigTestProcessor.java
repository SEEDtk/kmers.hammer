/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;

import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.java.erdb.DbConnection;
import org.theseed.proteins.hammer.HammerDb;

/**
 * This command analyzes a single FASTA file to look for hammer hits.  For each hit, it outputs the contig ID and location,
 * the ID of the feature from which the hammer was harvested, and the contig comment.  This information can be used to
 * determine how good the hits were.
 *
 *
 * The positional parameter is the name of the FASTA file.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 *
 * --hType		type of hammer database
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type
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

    // COMMAND-LINE OPTIONS

    /** batch size for database or web queries */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for queries")
    private int batchSize;

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
    @Option(name = "--parms", metaVar="user=xxx&pass=YYY", usage = "database parameter string (for MySQL")
    private String dbParms;

    /** input FASTA file name */
    @Argument(index = 0, metaVar = "input.fa", usage = "name of the FASTA file containing the sequences to search")
    private File inFile;

    @Override
    protected void setReporterDefaults() {
        this.batchSize = 200;
        this.hammerType = HammerDb.Type.MEMORY;
        this.dbEngine = DbConnection.Type.SQLITE;
        this.dbFile = null;
        this.dbUrl = null;
        this.dbParms = null;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1.");
        if (this.dbFile != null && ! this.dbFile.canRead())
            throw new FileNotFoundException("Database file " + this.dbFile + " is not found or unreadable.");
        // Set up the hammer database.  We must convert SQL exceptions for compatability.
        try {
            this.hammers = this.hammerType.create(this);
        } catch (SQLException e) {
            throw new IOException("Database error: " + e.toString());
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // TODO code for runReporter

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
