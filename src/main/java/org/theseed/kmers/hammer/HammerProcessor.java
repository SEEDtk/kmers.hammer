/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;
import java.util.TreeSet;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.hammer.GenomeHammerFactory;
import org.theseed.sequence.DnaKmers;
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
 * --para		if specified, parallel processing will be used wherever possible
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
    private TreeSet<String> genomeIdSet;
    /** role definition map for hammer roles */
    private RoleMap roleMap;

    // COMMAND-LINE OPTIONS

    /** kmer size */
    @Option(name = "--kmer", aliases = "-K", metaVar = "21", usage = "hammer DNA kmer size")
    private int kmerSize;

    /** optional filter file */
    @Option(name = "--filter", metaVar = "genomes.tbl", usage = "if specified, a file of genome IDs to select from the source")
    private File filterFile;

    /** TRUE for parallel mode */
    @Option(name = "--para", usage = "if specified, parallel processing will be used (default FALSE)")
    private boolean paraMode;

    /** genome source type */
    @Option(name = "--source", usage = "genome source type")
    private GenomeSource.Type sourceType;

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
        this.paraMode = false;
        this.sourceType = GenomeSource.Type.DIR;
        this.kmerSize = 20;
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
        // Set up parallelism and the kmer size.
        if (this.kmerSize <= 1)
            throw new ParseFailureException("Hammer kmer size must be at least 2.");
        DnaKmers.setKmerSize(this.kmerSize);
        GenomeHammerFactory.setParaMode(this.paraMode);
        if (this.paraMode)
            log.info("Parallel processing activated.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write the header line.
        writer.println("hammer\tfid");
        // Set up some counters.
        int gCount = 0;
        final int gTotal = this.genomeIdSet.size();
        // Loop through the genomes.
        long start = System.currentTimeMillis();
        for (String genomeId : this.genomeIdSet) {
            gCount++;
            // Compute the hammers.
            Genome genome = this.genomes.getGenome(genomeId);
            log.info("Processing genome {} of {}: {}.", gCount, gTotal, genome);
            GenomeHammerFactory factory = new GenomeHammerFactory(genome, this.roleMap);
            log.info("Enforcing precision rules.");
            factory.processFasta(this.contigFile);
            // Write out the hammers.
            factory.dumpHammers(writer);
            writer.flush();
            double seconds = (System.currentTimeMillis() - start) / (1000.0 * gCount);
            log.info("{} seconds/genome, estimated {} seconds remaining.", seconds, seconds * (gTotal - gCount));
        }
    }

}
