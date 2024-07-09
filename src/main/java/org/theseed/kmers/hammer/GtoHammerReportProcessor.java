/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.hammer.ScoreMap;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command will produce a report on the hammers in a set of genomes.  For each genome, we will count the
 * hammers using the specified method (COUNT or STRENGTH) and produce an output indicating the correlations between
 * each genome and the representative genomes indicated by the hammers.  The output report will have the ID and
 * name of each input genome and each representative genome as well as the hammer score.
 *
 * The positional parameters are the name of the repXX.stats.tbl file for the representative genomes used to
 * generate the hammers and the name (file or directory) of the input genome source.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 *
 * --source		type of genome source (default DIR)
 * --hType		type of hammer database (default MEMORY)
 * --method		voting method to use (default COUNT)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type
 * --para		use parallel processing to improve performance
 * --min		minimum score to report (default 10)
 *
 * @author Bruce Parrello
 *
 */
public class GtoHammerReportProcessor extends BaseHammerUsageProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GtoHammerReportProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** map of representative genome IDs to names */
    private Map<String, String> repGenMap;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of input genome source")
    private GenomeSource.Type sourceType;

    /** parallel-processing flag */
    @Option(name = "--para", usage = "if specified, parallel processing will be used")
    private boolean paraFlag;

    /** minimum score to report on output */
    @Option(name = "--min", metaVar = "2.0", usage = "minimum score to include on report")
    private double minScore;

    /** name of the representative-genome set stats file */
    @Argument(index = 0, metaVar = "repXX.stats.tbl", usage = "representative-genome stats file containing the repgen names",
            required = true)
    private File repStatsFile;

    /** name of the input genome source */
    @Argument(index = 1, metaVar = "inDir", usage = "input genome source file or directory", required = true)
    private File inDir;


    @Override
    protected void setHammerDefaults() {
        this.paraFlag = false;
        this.sourceType = GenomeSource.Type.DIR;
        this.minScore = 10.0;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Insure the genome source exists.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " is not found.");
        // Validate and load the rep-stats file.
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Repgen stats file " + this.repStatsFile + "is not found or unreadable.");
        this.repGenMap = TabbedLineReader.readMap(this.repStatsFile, "rep_id", "rep_name");
        log.info("{} representative genomes found in {}.", this.repGenMap.size(), this.repStatsFile);
        // Connect to the genome source.
        this.genomes = this.sourceType.create(this.inDir);
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        // Write the report header.
        writer.println("genome_id\tgenome_name\trep_id\trep_name\tscore\troles");
        // Get all the genome IDs for this source and build a stream.
        var genomeIDs = this.genomes.getIDs();
        Stream<String> genomeStream = this.makePara(genomeIDs.stream(), this.paraFlag);
        log.info("{} genomes to process from input source.", genomeIDs.size());
        // Process all the genomes.
        genomeStream.forEach(x -> this.processGenome(hammerDb, writer, x));
    }

    /**
     * Count the hammer hits in a single genome.
     *
     * @param hammerDb		hammer database to use
     * @param writer		output writer for the report
     * @param genomeID		ID of the genome to process
     */
    private void processGenome(HammerDb hammerDb, PrintWriter writer, String genomeID) {
        // Load the genome into memory.
        Genome genome = this.genomes.getGenome(genomeID);
        log.info("Processing genome {}.", genome);
        // We need to convert the contigs into a collection of sequences.
        Collection<Sequence> contigs = genome.getContigs().stream().map(x -> new Sequence(x.getId(), "", x.getSequence()))
                .collect(Collectors.toList());
        // Compute the scores.
        ScoreMap scoreMap = hammerDb.findClosest(contigs);
        // Write the results.
        this.writeResults(writer, genome, scoreMap);
    }

    /**
     * Write the hammer scores to the output report.
     *
     * @param writer		output report writer
     * @param genome		genome being scored
     * @param scoreMap		repgen scores from the hammers
     */
    private synchronized void writeResults(PrintWriter writer, Genome genome, ScoreMap scoreMap) {
        // Get the genome ID and name.
        String prefix = genome.getId() + "\t" + genome.getName();
        // Loop through the scores.
        int found = 0;
        List<ScoreMap.Count> scores = scoreMap.sortedCounts();
        for (var score : scores) {
            if (score.getCount() >= this.minScore) {
                writeCount(writer, prefix, score);
                found++;
            }
        }
        if (found <= 0) {
            // Nothing qualified.  If there was at least one score, show the first one.
            if (scores.size() > 0)
                this.writeCount(writer, prefix, scores.get(0));
            else {
                // Nothing hit.  Show a blank line.
                writer.println(prefix + "\t\t\t\t");
            }
        }
        writer.flush();
    }

    /**
     * Output a single report line.
     *
     * @param writer		output writer for report
     * @param prefix		prefix for the output line
     * @param score			score to output
     */
    private void writeCount(PrintWriter writer, String prefix, ScoreMap.Count score) {
        String repId = score.getKey();
        String repName = this.repGenMap.get(repId);
        writer.println(prefix + "\t" + repId + "\t" + repName + "\t" + Double.toString(score.getCount())
                + "\t" + Integer.toString(score.getNumRoles()));
    }

}
