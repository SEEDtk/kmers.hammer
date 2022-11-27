/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.hammer.GenomeHammerFactory;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.SequenceManager;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command will produce a table of hammers for all the genomes in a genome source.  A hammer is
 * a DNA kmer that is discriminating for key genes in a genome.  The hammer is considered precise if
 * it is not found in any other genomes of a representative set.  It is considered worthy if it
 * is found in genomes in the neighborhood of the original.  This command produces precise hammers but
 * does not test worthiness.
 *
 * This command requires a role definition file for the key genes (produced by HammerDefineProcessor),
 * a contig FASTA for all the genomes in the source (produced by by HammerFastaProcessor), and the
 * genome source for the genomes to be processed.  The genomes need to be at full detail level.
 * This is a ruinously slow process (a contig FASTA upwards of 32 gigabytes is not unheard of), so
 * it can be parallelized internally.  The input genomes can also be split up and processed in groups;
 * the "--filter" option is provided to facilitate this.
 *
 * The output file will be a two-column tab-delimited file with headers, the first column being the
 * hammer itself and the second being the ID of the feature that contained the hammer.  The default
 * output is to STDOUT.
 *
 * The positional parameters are the name (file or directory) of the genome source, the name of the
 * contig FASTA file, and the name of the role definition file.
 *
 * This is a long, slow command.  If an error or interruption occurs, you can resume processing by
 * specifying the output of a previous run using the "--resume" option.  In that case, the
 * previous file's hammers will be copied, and the genomes that were output will be skipped.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more detailed progress messages
 * -o	output file (if not STDOUT)
 * -K	hammer size (default 20)
 *
 * --source		genome source type (master, core, PATRIC, etc.); default is DIR
 * --filter		name of a filter file; if specified, must be a tab-delimited file with headers containing
 * 				genome IDs in the first column, and only the specified genomes will be processed
 * --para		if specified, parallel processing will be used to improve throughput
 * --resume		if specified, the name of the output file from a previous, interrupted run
 * --method		method for accessing the contigs (FILE to read during processing, MEM to keep in memory, default FILE)
 *
 * @author Bruce Parrello
 *
 */
public class HammerProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** filtered set of IDs for genomes to process */
    private Set<String> genomeIdSet;
    /** role definition map for hammer roles */
    private RoleMap roleMap;
    /** number of genomes processed */
    private int gCount;
    /** number of hammers found */
    private long hCount;
    /** sequence manager for contigs */
    private SequenceManager seqManager;
    /** start time of first genome */
    private long start;
    /** number of new genomes processed since the resume */
    private int newGenomes;
    /** total number of genomes in original set */
    private int nGenomes;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name = "--kmer", aliases = "-K", metaVar = "21", usage = "hammer DNA kmer size")
    private int kmerSize;

    /** optional filter file */
    @Option(name = "--filter", metaVar = "genomes.tbl", usage = "if specified, a file of genome IDs to select from the source")
    private File filterFile;

    /** TRUE for parallel mode */
    @Option(name = "--para", usage = "if specified, parallel processing will be used")
    private boolean paraMode;

    /** genome source type */
    @Option(name = "--source", usage = "genome source type")
    private GenomeSource.Type sourceType;

    /** sequence iteration method type */
    @Option(name = "--method", usage = "sequence iteration method")
    private SequenceManager.Type methodType;

    /** previous output file, for resume mode */
    @Option(name = "--resume", metaVar = "oldOutput.tbl", usage = "if specified, the name of a previous output file from an interrupted run to be resumed")
    private File resumeFile;

    /** name of the genome source file or directory */
    @Argument(index = 0, metaVar = "inDir", usage = "name of the genome source directory or file", required = true)
    private File inDir;

    /** name of the contig FASTA file */
    @Argument(index = 1, metaVar = "contigs.fna", usage = "name of the contig FASTA file for testing precision", required = true)
    private File contigFile;

    /** name of the role definition file */
    @Argument(index = 2, metaVar = "roles.for.hammers", usage = "name of the definition file for the hammer roles", required = true)
    private File roleFile;

    @Override
    protected void setReporterDefaults() {
        this.filterFile = null;
        this.resumeFile = null;
        this.paraMode = false;
        this.sourceType = GenomeSource.Type.DIR;
        this.kmerSize = 20;
        this.methodType = SequenceManager.Type.FILE;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure we can read the contig file.
        if (! this.contigFile.canRead())
            throw new FileNotFoundException("Contig FASTA file " + this.contigFile + " is not found or unreadable.");
        // Read in the role map.
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or unreadable.");
        this.roleMap = RoleMap.load(this.roleFile);
        // Set up the genome source.  First we attach the source itself.
        if (! this.inDir.exists())
            throw new FileNotFoundException("Genome source " + this.inDir + " not found.");
        log.info("Loading {} genome source {}.", this.sourceType, this.inDir);
        this.genomes = this.sourceType.create(this.inDir);
        // Handle the filter file.
        if (this.filterFile == null) {
            log.info("All {} genomes will be processed.", this.genomes.size());
            this.genomeIdSet = new TreeSet<String>(this.genomes.getIDs());
        } else if (! this.filterFile.canRead())
            throw new FileNotFoundException("Filter file " + this.filterFile + " is not found or unreadable.");
        else {
            Set<String> filter = TabbedLineReader.readSet(this.filterFile, "1");
            this.genomeIdSet = new TreeSet<String>();
            this.genomes.getIDs().stream().filter(x -> filter.contains(x)).forEach(x -> this.genomeIdSet.add(x));
            log.info("{} genomes remaining after filter using {}.", this.genomeIdSet.size(), this.filterFile);
        }
        // Check for the resume file.
        if (this.resumeFile != null && ! this.resumeFile.canRead())
            throw new FileNotFoundException("Previous-run output file " + this.resumeFile + " is not found or unreadable.");
        // Set up the kmer size.
        if (this.kmerSize <= 1)
            throw new ParseFailureException("Hammer kmer size must be at least 2.");
        DnaKmers.setKmerSize(this.kmerSize);
        // We don't parallelize at the low level, but on a genome basis instead.
        GenomeHammerFactory.setParaMode(false);
        // Set up the sequence manager.
        this.seqManager = this.methodType.create(this.contigFile);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Clear the genome and hammer counters.
        this.gCount = 0;
        this.hCount = 0;
        this.newGenomes = 0;
        this.nGenomes = this.genomeIdSet.size();
        // Write the header line.
        writer.println("hammer\tfid");
        // Handle the resume file, if any.
        if (this.resumeFile != null)
            this.processResume(writer);
        // Now form the genome IDs into an immutable list.  This gets us more efficient parallel processing.
        var genomeIdList = Collections.unmodifiableList(new ArrayList<String>(this.genomeIdSet));
        // Loop through the genomes.
        if (this.paraMode) {
            log.info("Parallel processing activated.");
            this.start = System.currentTimeMillis();
            genomeIdList.parallelStream().forEach(x -> this.processGenome(x, writer));
        } else {
            log.info("Serial processing used.");
            this.start = System.currentTimeMillis();
            genomeIdList.stream().forEach(x -> this.processGenome(x, writer));
        }
    }

    /**
     * Copy the previous run's output to the new output file, and remove the genomes already processed
     * from the genome ID set.
     *
     * @param writer	print stream for new output
     *
     * @throws IOException
     */
    private void processResume(PrintWriter writer) throws IOException {
        Set<String> removeSet = new HashSet<String>(1000);
        try (TabbedLineReader oldStream = new TabbedLineReader(this.resumeFile)) {
            log.info("Reading old output file {}.", this.resumeFile);
            for (var line : oldStream) {
                String hammer = line.get(0);
                String fid = line.get(1);
                String genome = Feature.genomeOf(fid);
                removeSet.add(genome);
                writer.println(hammer + "\t" + fid);
                this.hCount++;
                if (log.isInfoEnabled() && this.hCount % 100000 == 0)
                    log.info("{} hammers read.", this.hCount);
            }
        }
        log.info("{} genomes already processed with {} hammers.", removeSet.size(), this.hCount);
        this.genomeIdSet.removeAll(removeSet);
        log.info("{} genomes remaining to process.", this.genomeIdSet.size());
        this.gCount = removeSet.size();
    }

    private void processGenome(String genomeId, PrintWriter writer) {
        // Compute the hammers.
        Genome genome = this.genomes.getGenome(genomeId);
        log.info("Processing genome {}.", genome);
        GenomeHammerFactory factory = new GenomeHammerFactory(genome, this.roleMap);
        log.info("Enforcing precision rules.");
        int newHammers = 0;
        try {
            newHammers = factory.processFasta(this.seqManager);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        log.info("{} hammers found.", newHammers);
        // Write out the hammers.
        synchronized (this) {
            factory.dumpHammers(writer);
            writer.flush();
            this.gCount++;
            this.hCount += newHammers;
            if (log.isInfoEnabled()) {
                this.newGenomes++;
                Duration genomeTime = Duration.ofMillis(((System.currentTimeMillis() - this.start) / this.newGenomes + 1) * (this.nGenomes - this.gCount));
                log.info("{} genomes remaining, {} left.", this.nGenomes - this.gCount, genomeTime.toString());
            }
        }
    }

}
