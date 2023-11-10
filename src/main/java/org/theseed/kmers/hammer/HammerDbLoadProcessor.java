/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.erdb.utils.BaseDbProcessor;
import org.theseed.java.erdb.DbConnection;
import org.theseed.proteins.hammer.SqlHammerDb;

/**
 * This is a simple command that loads an ERDB-style hammer database from the flat file produced
 * using "HammerProcessor".
 *
 * The positional parameter is the name of the hammer input file.  The file must be tab-delimited, with
 * headers, containing the hammer sequence in the first column and the associated feature ID in the second.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 *
 * --type		type of database (default SQLITE)
 * --dbfile		database file name (SQLITE only)
 * --url		URL of database (host and name)
 * --parms		database connection parameter string (currently only MySQL)
 * --init		if specified, the name of an SQL file containing the database definition statements
 *
 * @author Bruce Parrello
 *
 */
public class HammerDbLoadProcessor extends BaseDbProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerDbLoadProcessor.class);

    // COMMAND-LINE OPTIONS

    /** SQL file to re-initialize the database */
    @Option(name = "--init", metaVar = "hammerDb.sql", usage = "if specified, an SQL file to redefine the database")
    private File initFile;

    /** input hammer file */
    @Argument(index = 0, metaVar = "hammers.tbl", usage = "hammer input file")
    private File inFile;

    @Override
    protected void setDbDefaults() {
        this.inFile = null;
        this.initFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Validate the initialization file.
        if (this.initFile != null && ! this.initFile.canRead())
            throw new FileNotFoundException("DB initialization SQL file " + this.initFile + " is not found or unreadable.");
        // Validate the input file.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Hammer input file " + this.inFile + " is not found or unreadable.");
        else
            log.info("Database will be loaded from {}.", this.inFile);
        return true;
    }

    @Override
    protected void runDbCommand(DbConnection db) throws Exception {
        if (this.initFile != null) {
            // Here we must initialize (or redefine) the database.
            log.info("Dropping current database tables.");
            db.clearTables();
            log.info("Initializing database using {}.", this.initFile);
            db.scriptUpdate(this.initFile);
            log.info("Database initialized.");
       }
        // Create the hammer database handler.  We set the query batch size to 100,
        // but it doesn't matter, since we are not doing any queries.
        SqlHammerDb hammerDb = new SqlHammerDb(db, 100);
        // Load the hammers.
        log.info("Loading from {}.", this.inFile);
        hammerDb.load(this.inFile);
    }

}
