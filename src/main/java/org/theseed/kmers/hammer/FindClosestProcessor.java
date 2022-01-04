/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Collection;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceDirectory;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;
import org.theseed.java.erdb.DbConnection;

/**
 * This method processes all the genomes in a directory and returns the closest genomes using
 * a hammer database.  The hammer database maps DNA kmers to feature IDs.  The feature ID can be
 * used to determine the source representative genome.
 *
 * The input directory can contain GTO files or FASTA files.  For a GTO file, the genome name and ID
 * will be pulled directly from the GTO.  For a FASTA file, the genome ID will be computed by looking
 * for a genome ID in the first header and the genome name will be the file name.
 *
 * There are two phases to this process. The hammer database is initialized in the first phase.  For
 * in-memory databases, this can be very slow, but the increased speed of the second phase makes up for
 * it.  Then the input directory is processed to create a list of sequences for each incoming GTO or
 * FASTA file.  The sequences are processed in parallel and the results written sequentially.
 *
 * The positional parameter is the name of the input directory.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -t	hammer type to be used in output headers
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
public class FindClosestProcessor extends BaseReportProcessor implements HammerDb.IParms, DbConnection.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FindClosestProcessor.class);
    /** hammer database */
    private HammerDb hammers;
    /** input directory streamer */
    private SequenceDirectory members;
    /** number of input members processed */
    private int memberCount;
    /** number of close-genomes not found */
    private int failureCount;

    // COMMAND-LINE OPTIONS

    /** batch size for database or web queries */
    @Option(name = "--batch", aliases = { "-b" }, metaVar = "100", usage = "batch size for queries")
    private int batchSize;

    /** type of hammer database */
    @Option(name = "--hType", usage = "type of hammer database")
    private HammerDb.Type hammerType;

    /** hammer title string */
    @Option(name = "--title", aliases = { "-t" }, metaVar = "hammer100",
            usage = "hammer type to appear in output file headers")
    private String hammerTitle;

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

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing sequence files", required = true)
    private File inDir;

    @Override
    protected void setReporterDefaults() {
        this.hammerTitle = "hammer";
        this.dbEngine = DbConnection.Type.SQLITE;
        this.dbFile = null;
        this.dbUrl = null;
        this.dbParms = null;
        this.batchSize = 4000;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Verify that the sequence directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.members = new SequenceDirectory(this.inDir);
        log.info("{} members found in {}.", this.members.size(), this.inDir);
        // Validate the batch size.
        if (this.batchSize < 1)
            throw new ParseFailureException("Batch size must be at least 1.");
        // We must convert SQL exceptions for compatability.Mo
        try {
            this.hammers = this.hammerType.create(this);
        } catch (SQLException e) {
            throw new IOException("Database error: " + e.toString());
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the output heading.
        writer.format("genome_id\tgenome_name\t%s.DNA.closest_genome1\t%s.DNA.closeness1\t%s.DNA.closest_genome2\t%s.DNA.closeness2%n",
                this.hammerTitle, this.hammerTitle, this.hammerTitle, this.hammerTitle);
        // Initialize the member counter.
        this.memberCount = 0;
        // Process all the members.
        Stream<SequenceDirectory.Member> inStream = this.members.stream();
        long start = System.currentTimeMillis();
        inStream.forEach(x -> this.analyzeMember(x, writer));
        if (log.isInfoEnabled() && this.memberCount > 0) {
            double rate = (System.currentTimeMillis() - start) / (1000.0 * this.memberCount);
            log.info("{} members processed, at {} seconds per member. {} failures.", this.memberCount, rate,
                    this.failureCount);
        }
    }

    /**
     * Analyze the sequences of a member and output the closest genomes.
     *
     * @param member	member to analyze
     * @param writer	print writer for output
     */
    private void analyzeMember(SequenceDirectory.Member member, PrintWriter writer) {
        log.info("Processing {}: {}.", member.getId(), member.getName());
        Collection<Sequence> seqs = member.getSequences();
        CountMap<String> closeMap = this.hammers.findClosest(seqs);
        // We will build our output line in here.
        StringBuffer buffer = new StringBuffer(200);
        buffer.append(member.getId()).append("\t").append(member.getName());
        var closestList = closeMap.sortedCounts();
        // Insure we only show the two closest genomes.
        if (closestList.size() > 2)
            closestList = closestList.subList(0, 2);
        // Display the closest genomes found.
        for (CountMap<String>.Count counter : closestList)
            buffer.append("\t").append(counter.getKey()).append("\t").append(counter.getCount());
        if (closestList.size() < 2) {
            // Here we have a failure condition.
            this.recordFailure();
            // Fill out any empty columns.
            for (int i = closestList.size(); i < 2; i++)
                buffer.append("\t\t");
        }
        // Now we write this synchronized.
        this.write(writer, buffer.toString());
    }

    /**
     * Denote that we failed to find two close genomes.
     */
    private synchronized void recordFailure() {
        this.failureCount++;
    }

    /**
     * Write a line of output for a member.  This method is synchronized, to prevent interference between
     * threads.
     *
     * @param writer		output print writer
     * @param string		data line for the current member
     */
    private synchronized void write(PrintWriter writer, String string) {
        writer.println(string);
        this.memberCount++;
        log.info("{} of {} members processed. {} failures.", this.memberCount, this.members.size(), this.failureCount);
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
