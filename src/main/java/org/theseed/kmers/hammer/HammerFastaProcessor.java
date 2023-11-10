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
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.proteins.kmers.reps.RepGenomeDb;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This command reads all the genomes from a genome source and converts them into a single DNA FASTA file containing
 * the contigs.  The label of each contig will be based on the contig ID, and the comment will be
 * the genome name, the ID of the closest representative, and the distance.
 *
 * The positional parameters are the name of the genome source directory (or file, in the case of a PATRIC source),
 * the name of the representative-genome database, and the name of the output file.
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
    /** representative-genome database */
    private RepGenomeDb repDb;

    // COMMAND-LINE OPTIONS

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** input file or directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input genome directory or file", required = true)
    private File inDir;

    /** representative-genome database */
    @Argument(index = 1, metaVar = "repXXX.ser", usage = "representative-genome database", required = true)
    private File repDbFile;

    /** output file name */
    @Argument(index = 2, metaVar = "outFile.fna", usage = "output FASTA file name", required = true)
    private File outFile;

    @Override
    protected void setDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Load the representative-genome database.
        if (! this.repDbFile.canRead())
            throw new FileNotFoundException("Repgen database file " + this.repDbFile + " is not found or unreadable.");
        log.info("Loading repgen database from {}.", this.repDbFile);
        this.repDb = RepGenomeDb.load(this.repDbFile);
        // Set up the input genome source.
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
                    // Here we need to compute the comment.  To do that, we need to find the closest
                    // representative.  In rare cases, there will be none, in which case we use an empty string.
                    var seedProt = this.repDb.getSeedProtein(genome);
                    var representation = this.repDb.findClosest(seedProt);
                    String closestId = representation.getGenomeId();
                    if (closestId == null) closestId = "";
                    var comment = genome.getName() + "\t" + closestId + "\t" + representation.getDistance();
                    // Form the output contig.
                    Sequence seq = new Sequence(label, comment, contig.getSequence());
                    outStream.write(seq);
                    contigCount++;
                }
            }
            log.info("{} contigs written for {} genomes.", contigCount, count);
            outStream.flush();
        }
    }

}
