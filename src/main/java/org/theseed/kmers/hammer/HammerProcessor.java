/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeUniSequences;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command finds all the hammers for a specified set of genomes.  A hammer is defined in this
 * case as a DNA kmer that occurs in one genome of the set but not any of the others.  We restrict our
 * search to a particular set of universal proteins found in most genomes.  This is done to keep the
 * hammer set to a manageable size.
 *
 * The output file will be a two-column table consisting of the hammer sequence followed by the ID of
 * the feature than contained it. The file will have a header line.
 *
 * The program operates in two phases.  In the first phase, we will load the DNA sequences of the
 * universal proteins into memory as sequence sets.  The second phase iterates through these sequence
 * sets in parallel.  For each set, it builds a kmer hash and then removes any kmers found in the
 * other sequences.  The level of parallelism needs to be managed to prevent memory overload:  the
 * kmer hash is significantly larger than the set of sequences.
 *
 * Once it reaches the second phase, the program is restartable.  If we restart, we track the genomes
 * previously processed and skip those when we rerun phase two.
 *
 * The positional parameters are the name of the role definition file, the directory or file of the
 * input genome source, and the name of the output file.  If the output file exists, we will assume
 * that we are resuming.  Use the "--clear" option to override this behavior.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -K	kmer size (default 20)
 *
 * --clear	erase the output file before starting
 * --type	type of genome source (default DIR)
 *
 * @author Bruce Parrello
 *
 */
public class HammerProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerProcessor.class);
    /** map of genome IDs to universal-role descriptors */
    private Map<String, GenomeUniSequences> genomeSequences;
    /** number of hammers found */
    private long hammerCount;
    /** number of genomes completed */
    private int genomeCount;

    // COMMAND-LINE OPTIONS

    /** dna kmer size */
    @Option(name = "--kmer", aliases = { "-K" }, metaVar = "15", usage = "dna kmer size")
    private int kmerSize;

    /** TRUE to clear the output file before starting */
    @Option(name = "--clear", usage = "if specified, the output file will be erased before starting")
    private boolean clearFlag;

    /** type of genome source (PATRIC, DIR, MASTER) */
    @Option(name = "--type", usage = "type of genome source (PATRIC, DIR, MASTER)")
    private GenomeSource.Type sourceType;

    /** role definition file name */
    @Argument(index = 0, metaVar = "role.definitions", usage = "name of file containing the definitions for the roles of interest")
    private File roleFile;

    /** genome input source */
    @Argument(index = 1, metaVar = "genomeDir", usage = "name of file or directory containing source genomes")
    private File genomeDir;

    /** output file */
    @Argument(index = 2, metaVar = "outFile", usage = "output file name")
    private File outFile;

    @Override
    protected void setDefaults() {
        this.kmerSize = 20;
        this.clearFlag = false;
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the kmer size is valid.
        if (this.kmerSize < 4)
            throw new ParseFailureException("Kmer size must be at least 4.");
        GenomeUniSequences.setKmerSize(this.kmerSize);
        // Validate the input files.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("Input genome source " + this.genomeDir + " does not exist.");
        if (! this.roleFile.canRead())
            throw new FileNotFoundException("Role definition file " + this.roleFile + " is not found or is unreadable.");
        // If the clear flag is set, delete the output file.
        if (this.outFile.exists()) {
            if (this.clearFlag) {
                log.info("Deleting old output file.");
                FileUtils.forceDelete(this.outFile);
            } else if (! this.outFile.canRead())
                throw new FileNotFoundException("Old output file " + this.outFile + " is unreadable.");
            else if (! this.outFile.canWrite())
                throw new FileNotFoundException("Output file " + this.outFile + " is unwritable.");
        }
        // Denote no hammers have been found yet.
        this.hammerCount = 0;
        this.genomeCount = 0;
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // We do most of the work in subroutines so that memory is managed more efficiently.
        // First, we build the genome sequences.  This part uses the genome source and the role map.
        long start = System.currentTimeMillis();
        this.genomeSequences = this.buildGenomeDescriptorMap();
        if (log.isInfoEnabled()) {
            Duration length = Duration.ofMillis(System.currentTimeMillis() - start);
            log.info("{} to build genome sequence map.", length.toString());
        }
        // Now, we check the output file to build a list of the descriptors to process.  This part
        // uses the set of genomes already in the output file.
        Collection<GenomeUniSequences> genomes = this.selectGenomes();
        log.info("{} genomes remaining to process.", genomes.size());
        // Finally, we open the output file and process the genomes in parallel.  A synchronized
        // subroutine is used to write the output.
        start = System.currentTimeMillis();
        int startingGenomes = this.genomeCount;
        try (PrintWriter writer = this.getOutputStream()) {
            genomes.parallelStream().forEach(x -> this.findHammers(x, writer));
        }
        if (log.isInfoEnabled()) {
            double rate = (System.currentTimeMillis() - start) / (1000.0 * (this.genomeCount - startingGenomes));
            log.info("{} total hammers found in {} genomes at {} seconds/genome", this.hammerCount, this.genomeCount, rate);
        }
    }

    /**
     * @return a map from genome IDs to universal-protein DNA sequence maps
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private Map<String, GenomeUniSequences> buildGenomeDescriptorMap() throws IOException, ParseFailureException {
        // Read in the role map.
        RoleMap roleMap = RoleMap.load(this.roleFile);
        // Get the genome source.
        log.info("Scanning genome source of type {} in {}.", this.sourceType, this.genomeDir);
        GenomeSource genomes = this.sourceType.create(this.genomeDir);
        int n = genomes.size();
        log.info("{} genomes found.", n);
        var retVal = new HashMap<String, GenomeUniSequences>((n + 1) * 4 / 3);
        // Loop through the genomes, creating the sequence objects.
        int gCount = 1;
        for (Genome genome : genomes) {
            log.info("Scanning genome #{} of {}: {}.", gCount, n, genome);
            retVal.put(genome.getId(), new GenomeUniSequences(genome, roleMap));
            gCount++;
        }
        return retVal;
    }

    /**
     * @return the universal-protein DNA sequence maps that still need to be processed
     *
     * @throws IOException
     */
    private Collection<GenomeUniSequences> selectGenomes() throws IOException {
        Collection<GenomeUniSequences> retVal = null;
        // Do we have an old output file?
        if (! this.outFile.exists()) {
            // No.  Return the whole set.
            retVal = this.genomeSequences.values();
            log.info("New output file.  All genomes will be processed.");
            this.genomeCount = 0;
        } else {
            // Here we must scan the old output file for already-processed genomes.
            Set<String> oldGenomes = new TreeSet<String>();
            log.info("Scanning old output file for progress.");
            try (TabbedLineReader oldStream = new TabbedLineReader(this.outFile)) {
                for (TabbedLineReader.Line line : oldStream) {
                    // The second column is the feature ID. We extract its genome ID.
                    String genome = Feature.genomeOf(line.get(1));
                    oldGenomes.add(genome);
                    // Count the hammer.
                    this.hammerCount++;
                }
            }
            log.info("{} genomes already found in output file.", oldGenomes.size());
            this.genomeCount = oldGenomes.size();
            // Form a set of universal-protein DNA sequence maps exclusing the old genomes.
            retVal = this.genomeSequences.values().stream()
                    .filter(x -> ! oldGenomes.contains(x.getGenomeId()))
                    .collect(Collectors.toList());
        }
        return retVal;
    }

    /**
     * This method determines whether we are creating the output file or appending to it,
     * and opens it accordingly.
     *
     * @return the print writer for the output stream
     * @throws IOException
     */
    private PrintWriter getOutputStream() throws IOException {
        PrintWriter retVal;
        if (this.outFile.exists()) {
            log.info("Appending new results to {}.", this.outFile);
            FileWriter appender = new FileWriter(this.outFile, true);
            retVal = new PrintWriter(appender);
        } else {
            log.info("Writing to new file {}.", this.outFile);
            retVal = new PrintWriter(this.outFile);
            retVal.println("kmer\tfid");
            retVal.flush();
        }
        return retVal;
    }

    /**
     * This method finds the hammers in the specified genome, and writes them to the output
     * stream.
     *
     * @param genome	universal-protein DNA sequence map for the genome to process
     * @param writer 	output writer on which to list the hammers
     */
    private void findHammers(GenomeUniSequences genome, PrintWriter writer) {
        // Get the kmer map for this genome.
        Map<String, String> kmerMap = genome.getKmerMap();
        log.info("{} kmers found in {}.", kmerMap.size(), genome);
        // Loop through the other genomes, updating the kmer map.
        for (GenomeUniSequences other : this.genomeSequences.values()) {
            // Only compare against the other genomes.
            if (! genome.equals(other))
                other.processKmers(kmerMap);
        }
        // Output the hammers.
        int found = this.writeHammers(kmerMap, writer);
        log.info("{} hammers found in {}.", found, genome);
    }

    /**
     * This method will write the hammers found to the output stream.  It also
     * increments the hammer counter.
     *
     * @param kmerMap	map of hammer sequences to feature IDs
     * @param writer	output stream to receive the hammers
     */
    private synchronized int writeHammers(Map<String, String> kmerMap, PrintWriter writer) {
        // Update the hammer count.  It is safe to do here because we are synchronized.  Note
        // that there are two real hammers for every hammer sequence.
        int retVal = kmerMap.size() * 2;
        this.hammerCount += retVal;
        this.genomeCount++;
        // Write the hammers.  We flush at the end because we want all the hammers for a genome
        // to be output together, insuring the restart works properly.
        for (Map.Entry<String, String> kmerEntry : kmerMap.entrySet()) {
            String fid = kmerEntry.getValue();
            String seq = kmerEntry.getKey();
            writer.println(seq + "\t" + fid);
            writer.println(Contig.reverse(seq) + "\t" + fid);
        }
        writer.flush();
        log.info("{} genomes processed.", this.genomeCount);
        return retVal;
    }

}
