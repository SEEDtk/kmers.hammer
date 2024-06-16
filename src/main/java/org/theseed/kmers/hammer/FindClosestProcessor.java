/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.hammer.ScoreMap;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceDirectory;
import org.theseed.utils.BaseHammerUsageProcessor;

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
 * --hType		type of hammer database (default MEMORY)
 * --method		voting method to use (default COUNT)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type
 * --para		process the sequences in parallel
 *
 * @author Bruce Parrello
 *
 */
public class FindClosestProcessor extends BaseHammerUsageProcessor {

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
    /** process time spent computing closest-genomes */
    private long processTime;

    // COMMAND-LINE OPTIONS

    /** hammer title string */
    @Option(name = "--title", aliases = { "-t" }, metaVar = "hammer100",
            usage = "hammer type to appear in output file headers")
    private String hammerTitle;

    /** TRUE to turn on parallelization */
    @Option(name = "--para", usage = "if specified, multiple genomes will be run at the same time")
    private boolean paraFlag;

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory containing sequence files", required = true)
    private File inDir;

    @Override
    protected void setHammerDefaults() {
        this.hammerTitle = "hammer";
        this.paraFlag = false;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Verify that the sequence directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        this.members = new SequenceDirectory(this.inDir);
        log.info("{} members found in {}.", this.members.size(), this.inDir);
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        // Save the hammer database.
        this.hammers = hammerDb;
        // Write the output heading.
        writer.format("genome_id\tgenome_name\t%s.DNA.closest_genome1\t%s.DNA.closeness1\t%s.DNA.closest_genome2\t%s.DNA.closeness2%n",
                this.hammerTitle, this.hammerTitle, this.hammerTitle, this.hammerTitle);
        // Initialize the member counter and process time.
        this.memberCount = 0;
        this.processTime = 0;
        // Process all the members.
        Stream<SequenceDirectory.Member> inStream = this.members.stream(this.paraFlag);
        long start = System.currentTimeMillis();
        inStream.forEach(x -> this.analyzeMember(x, writer));
        if (log.isInfoEnabled() && this.memberCount > 0) {
            double rate = (System.currentTimeMillis() - start) / (1000.0 * this.memberCount);
            log.info("{} members processed, at {} seconds per member. {} failures.", this.memberCount, rate,
                    this.failureCount);
            // Compute the actual scan time per member.
            double scanRate = this.processTime / (this.memberCount * 1000.0);
            log.info("{} microseconds per genome to scan.", scanRate);
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
        long memberTimer = System.nanoTime();
        ScoreMap closeMap = this.hammers.findClosest(seqs);
        memberTimer = System.nanoTime() - memberTimer;
        // We will build our output line in here.
        StringBuffer buffer = new StringBuffer(200);
        buffer.append(member.getId()).append("\t").append(member.getName());
        var closestList = closeMap.sortedCounts();
        // Insure we only show the two closest genomes.
        if (closestList.size() > 2)
            closestList = closestList.subList(0, 2);
        // Display the closest genomes found.
        for (ScoreMap.Count counter : closestList)
            buffer.append("\t").append(counter.getKey()).append("\t").append(counter.getCount());
        if (closestList.size() < 2) {
            // Here we have a failure condition.
            this.recordFailure();
            // Fill out any empty columns.
            for (int i = closestList.size(); i < 2; i++)
                buffer.append("\t\t");
        }
        // Now we write this synchronized.
        this.write(writer, buffer.toString(), memberTimer);
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
     * @param memberTimer	nanoseconds spent analyzing the member
     */
    private synchronized void write(PrintWriter writer, String string, long memberTimer) {
        writer.println(string);
        this.memberCount++;
        this.processTime += memberTimer;
        log.info("{} of {} members processed. {} failures.", this.memberCount, this.members.size(), this.failureCount);
    }

}
