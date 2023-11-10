/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command processes a FASTA files against a set of hammers.  Hammers that are found in the source
 * sequences will be removed from the set.  This is an expensive process.  First, we load the hammer database
 * and process the sequences from the source.  This gives us a set of hammers that hit.  Then we read the hammer
 * load file and filter out the hammers in the set.
 *
 * The positional parameter is the name of the genome source directory or file.  The new hammer load file will
 * be written to the standard output.
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
 * --loader		name of database load file (if database type is SQL)
 *
 * @author Bruce Parrello
 *
 */
public class HammerFilterProcessor extends BaseHammerUsageProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerFilterProcessor.class);
    /** set of hammers to filter out */
    private Set<String> badHammerSet;
    /** FASTA files to process */
    private File[] fastaFiles;
    /** name of hammer load file */
    private File hammerLoadFile;
    /** file filter for FASTA files */
    private FilenameFilter FASTA_FILTER = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return name.endsWith(".fa") || name.endsWith(".fna") || name.endsWith(".fasta");
        }
    };


    // COMMAND-LINE OPTIONS

    /** load file name (only needed for SQL databases) */
    @Option(name = "--loader", metaVar = "hammers.tbl", usage = "hammer load file name (for SQL databases only)")
    private File loadFile;

    /** input FASTA file directory */
    @Argument(index = 0, metaVar = "inDir", usage = "FASTA file directory")
    private File inDir;

    @Override
    protected void setHammerDefaults() {
        this.loadFile = null;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input genome source " + this.inDir + " is not found.");
        // Get the hammer load file.
        this.hammerLoadFile = this.getLoadFileName();
        if (this.hammerLoadFile == null)
            this.hammerLoadFile = this.loadFile;
        if (this.hammerLoadFile == null)
            throw new ParseFailureException("Load file name must be specified for this type of hammer database.");
        else if (! this.hammerLoadFile.canRead())
            throw new FileNotFoundException("Hammer load file " + this.hammerLoadFile + " is not found or unreadable.");
        else
            log.info("Hammer load file is {}.", this.hammerLoadFile);
        // Now find the FASTA files.
        this.fastaFiles = this.inDir.listFiles(FASTA_FILTER);
        log.info("{} FASTA files found in {}.", this.fastaFiles.length, this.inDir);
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        // Create the bad-hammer set.
        this.badHammerSet = ConcurrentHashMap.newKeySet();
        // Process each file individually.  We will do the sequences in parallel, however.
        for (File fastaFile : this.fastaFiles) {
            log.info("Processing file {}.", fastaFile);
            int oldSize = this.badHammerSet.size();
            // Load all the sequences.
            List<Sequence> seqs = FastaInputStream.readAll(fastaFile);
            log.info("{} contigs found in file.", seqs.size());
            seqs.parallelStream().forEach(x -> this.processSequence(hammerDb, x));
            log.info("{} new bad-hammers found in {}.", this.badHammerSet.size() - oldSize, fastaFile);
        }
        log.info("{} total bad hammers found.", this.badHammerSet.size());
        // Open the load file.
        try (TabbedLineReader loadStream = new TabbedLineReader(this.hammerLoadFile)) {
            // Write the header line.
            writer.println(loadStream.header());
            // Loop through the data lines.  The hammer itself is always in the first column.
            int count = 0;
            for (var line : loadStream) {
                String hammer = line.get(0);
                if (! this.badHammerSet.contains(hammer)) {
                    writer.println(line.toString());
                    count++;
                }
            }
            log.info("{} hammers output.", count);
        }
    }

    /**
     * Search for bad hammers in the specified sequence.  We search the sequence on both strands, and add
     * the bad hammers to the bad-hammer set.
     *
     * @param hammerDb	hammer database
     * @param sequence	sequence to search
     */
    private void processSequence(HammerDb hammerDb, Sequence sequence) {
        log.info("Processing sequence {} of length {}.", sequence.getLabel(), sequence.length());
        String dna = sequence.getSequence();
        this.checkSequence(hammerDb, dna);
        dna = Contig.reverse(dna);
        this.checkSequence(hammerDb, dna);
    }

    /**
     * Search for bad hammers in a single sequence.  This method does not bother with strands.
     *
     * @param hammerDb	hammer database
     * @param sequence	DNA sequence to check
     */
    private void checkSequence(HammerDb hammerDb, String sequence) {
        var hammersFound = hammerDb.findHammers(sequence);
        this.badHammerSet.addAll(hammersFound);
    }

}
