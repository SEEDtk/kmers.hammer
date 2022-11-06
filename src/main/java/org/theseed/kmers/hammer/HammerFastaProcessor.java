/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command reads all the genomes from a genome source and converts them into a single DNA FASTA file containing
 * the contigs.  The label of each contig will be based on the contig ID, and the comment will be
 * the genome ID itself.
 *
 * The positional parameters are the name of the genome source directory (or file, in the case of a PATRIC source), and
 * the name of the output file.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent error messages
 *
 * --source		type of genome source (master, PATRIC, directory); default DIR
 *
 * @author Bruce Parrello
 *
 */
public class HammerFastaProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerFastaProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** set of contig IDs already used */
    private Set<String> contigIds;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** input file or directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input genome directory or file", required = true)
    private File inDir;

    /** output file name */
    @Argument(index = 1, metaVar = "outFile.fna", usage = "output FASTA file name", required = true)
    private File outFile;

    @Override
    protected void setDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        if (! this.inDir.exists())
            throw new FileNotFoundException("Input genome source " + this.inDir + " not found.");
        log.info("Loading genome source from {}.", this.inDir);
        this.genomes = this.sourceType.create(this.inDir);
        // Create the contig ID set.  This is used to prevent duplicate labels.
        this.contigIds = new HashSet<String>(this.genomes.size() * 125);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Open the output stream.
        try (FastaOutputStream outStream = new FastaOutputStream(this.outFile)) {
            // Loop through the genomes.
            final int total = this.genomes.size();
            int count = 0;
            int contigCount = 0;
            for (Genome genome : this.genomes) {
                count++;
                log.info("Processing {} of {}: {}.", count, total, genome);
                String genomeId = genome.getId();
                // Loop through the contigs.
                for (Contig contig : genome.getContigs()) {
                    String label = contig.getId();
                    if (this.contigIds.contains(label)) {
                        // To insure the label is unique, suffix the genome ID.
                        label += "_" + genomeId;
                    } else
                        this.contigIds.add(label);
                    Sequence seq = new Sequence(label, genomeId, contig.getSequence());
                    outStream.write(seq);
                    contigCount++;
                }
            }
            log.info("{} contigs written for {} genomes.", contigCount, count);
            outStream.flush();
        }
    }

}
