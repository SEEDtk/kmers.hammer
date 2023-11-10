/**
 *
 */
package org.theseed.utils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.basic.ParseFailureException;
import org.theseed.java.erdb.DbConnection;
import org.theseed.kmers.hammer.FindClosestProcessor;

/**
 * This is the base class for all reports produced using the hammers.
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
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseHammerUsageProcessor extends BaseReportProcessor implements HammerDb.IParms, DbConnection.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FindClosestProcessor.class);
    /** hammer database */
    private HammerDb hammers;

    // COMMAND-LINE OPTIONS

    /** batch size for database or web queries */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for queries")
    private int batchSize;

    /** type of hammer database */
    @Option(name = "--hType", usage = "type of hammer database")
    private HammerDb.Type hammerType;

    /** voting method to use */
    @Option(name = "--method", usage = "voting method for determining closest genome")
    private HammerDb.Method methodType;

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

    @Override
    final protected void setReporterDefaults() {
        this.dbEngine = DbConnection.Type.SQLITE;
        this.dbFile = null;
        this.dbUrl = null;
        this.dbParms = null;
        this.batchSize = 4000;
        this.hammerType = HammerDb.Type.MEMORY;
        this.methodType = HammerDb.Method.COUNT;
        this.setHammerDefaults();
    }

    /**
     * Specify the default values for subclass options.
     */
    protected abstract void setHammerDefaults();

    @Override
    final protected void validateReporterParms() throws IOException, ParseFailureException {
        // Validate the batch size.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1.");
        // We must convert SQL exceptions for compatability.
        try {
            this.hammers = this.hammerType.create(this);
            this.hammers.setMethod(this.methodType);
        } catch (SQLException e) {
            throw new IOException("Database error: " + e.toString());
        }
        // Validate the subclass parameters.
        this.validateHammerParms();
    }

    /**
     * Validate the subclass parameters.
     *
     * @throws IOException
     * @throws ParseFailureException
     */
    protected abstract void validateHammerParms() throws IOException, ParseFailureException;

    @Override
    final protected void runReporter(PrintWriter writer) throws Exception {
        this.runHammers(this.hammers, writer);
    }

    /**
     * Process the hammers against the input.
     *
     * @param hammerDb	hammer database
     * @param writer	output stream for report
     */
    protected abstract void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception;

    /**
     * @return TRUE if the hammer database knows its load file name
     */
    public boolean isLoadFileKnown() {
        return this.hammerType.isLoadFileKnown();
    }

    /**
     * @return the name of the load file, or NULL if it is unknown
     */
    public File getLoadFileName() {
        return this.hammers.getLoadFile();
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

    /**
     * @return the count method for this hammer operation
     */
    public HammerDb.Method getMethod() {
        return this.hammers.getCountMethod();
    }

}
