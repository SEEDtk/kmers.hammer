/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.hammer.SequenceKmerIterable;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command reads a bad-hits file from ReadTestProcessor and looks for hammers that are too close to sequences in
 * the expected genome or in the hammer database at large.  A list of these suspicious hammers is output.  If we find
 * a lot of them, we will consider adding pruning steps to the hammer database.
 *
 * A "close sequence" here is one that differs by a single nucleotide.  For any hammer, there are K*4 close sequences,
 * when K is the hammer size.  If a close sequence is found in the expected representative (which will be very different
 * from the genome the hammer was harvested from), then we assume the hammer is mutating fast enough to create a false
 * positive.  If the same sequence is also found in the hammer database, then we can argue that searching the hammer
 * database, rather than all the representative genomes, is sufficient to remove false positives of this sort.
 *
 * We do this in two passes.  First, we create a reduced file that has one instance of each hammer found along with
 * its total weight and the corresponding expected-genome set.  Then we process this reduced file against the expected-genomes
 * and the full hammer database.
 *
 * The positional parameters should be the name of the original report file containing the hammers to analyze, the
 * name of the blacklist file showing the expected genomes, and the name of a directory containing the repgen GTOs.
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
 * --gtoType	type of genome source for the representative genomes (default DIR)
 * --temp		name of a temporary file to contain the reduced bad-hits data (default: "reduced.temp.tbl" in the current directory)
 * --keep		if specified, the temporary file will not be deleted when the command completes
 *
 * @author Bruce Parrello
 *
 */
public class BadCheckProcessor extends BaseHammerUsageProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BadCheckProcessor.class);
    /** genome source for representative genomes */
    private GenomeSource repGenomes;
    /** map of sample IDs to expected genome sets */
    private Map<String, Set<String>> blackMap;
    /** repgen cache */
    private Map<String, Genome> repGtos;
    /** empty expected=genome set */
    private Set<String> NULL_SET = Collections.emptySet();
    /** array of nucleic acids */
    private char[] ACIDS = new char[] { 'a', 'c', 'g', 't' };

    // COMMAND-LINE OPTIONS

    /** type of repgen genome source */
    @Option(name = "--gtoType", usage = "repgen source type")
    private GenomeSource.Type gtoType;

    /** name of the temporary file for the reduced bad-hit data */
    @Option(name = "--temp", metaVar = "temp.tbl", usage = "name of the temporary file for the reduced hit data")
    private File tempFile;

    /** TRUE if we should preserve the temporary file */
    @Option(name = "--keep", usage = "If specified, the temporary file will not be deleted on exit")
    private boolean keepFlag;

    /** name of bad-hits input file */
    @Argument(index = 0, metaVar = "badHits.tbl", usage = "output file from readTest listing the hits to analyze", required = true)
    private File inFile;

    /** name of blacklist file */
    @Argument(index = 1, metaVar = "blacklist.tbl", usage = "file mapping sample IDs to expected genomes.", required = true)
    private File blackFile;

    /** directory or file name of repgen source */
    @Argument(index = 2, metaVar = "gtoDir", usage = "file name or directory of representative genome source", required = true)
    private File gtoDir;

    /**
     * This utility object describes a hammer of interest.  It contains the total weight and the set of expected genomes.
     */
    protected static class HammerData {

        /** total weight for this hammer */
        private double weight;
        /** set of expected genomes for this hammer */
        private Set<String> expected;

        /**
         * Create a blank, empty hammer-data descriptor.
         */
        public HammerData() {
            this.weight = 0.0;
            this.expected = new TreeSet<String>();
        }

        /**
         * Update this hammer-data descriptor.
         *
         * @param newWeight		weight to add
         * @param repSet		set of new repgen IDs to add
         */
        public void update(double newWeight, Set<String> repSet) {
            this.weight += newWeight;
            this.expected.addAll(repSet);
        }

        /**
         * @return the total weight for this hammer
         */
        public double getWeight() {
            return this.weight;
        }

        /**
         * @return the set of expected repgen IDs
         */
        public Set<String> getExpected() {
            return this.expected;
        }

    }

    @Override
    protected void setHammerDefaults() {
        this.keepFlag = false;
        this.tempFile = new File(System.getProperty("user.dir"), "reduced.temp.tbl");
        this.gtoType = GenomeSource.Type.DIR;
    }

    @Override
    protected void validateHammerParms() throws IOException, ParseFailureException {
        // Load the black list.
        this.blackMap = this.readBlackList(this.blackFile);
        // Insure the input file and the repgen source exist.
        if (! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or invalid.");
        if (! this.gtoDir.exists())
            throw new FileNotFoundException("Repgen source " + this.gtoDir + " does not exist.");
        // Check to insure all the blacklisted genomes are in the genome source.
        log.info("Verifying genome source {} (type {}) against the black list in {}.", this.gtoDir, this.gtoType, this.blackFile);
        this.repGenomes = this.gtoType.create(this.gtoDir);
        Set<String> repIDs = this.repGenomes.getIDs();
        for (Set<String> blackSet : this.blackMap.values()) {
            for (String repId : blackSet) {
                if (! repIDs.contains(repId))
                    throw new ParseFailureException("Genome source " + this.gtoDir + " is missing required genome " + repId + ".");
            }
        }
        // Now we need to load the input file and create the temporary file.  Note we delete the temporary
        // file if we crash before closing it.
        boolean closed = false;
        try (PrintWriter tempStream = new PrintWriter(this.tempFile);
                TabbedLineReader inStream = new TabbedLineReader(this.inFile)) {
            Map<String, HammerData> hammerMap = this.readInFile(inStream);
            this.writeTempFile(tempStream, hammerMap);
            closed = true;
        } finally {
            if (! closed && ! this.keepFlag) {
                log.info("Deleting temporary file.");
                FileUtils.forceDelete(this.tempFile);
            }
        }
        // Create the GTO hash.
        this.repGtos = new HashMap<String, Genome>(100);
    }

    /**
     * Create a map of hammer-data descriptors from the input file.
     *
     * @param inStream		input stream
     *
     * @return a map from hammer strings to hammer-data descriptors
     *
     * @throws IOException
     */
    private Map<String, HammerData> readInFile(TabbedLineReader inStream) throws IOException {
        Map<String, HammerData> retVal = new HashMap<String, HammerData>();
        // Get the sample-id, hammer, and weight column indices.
        int sampColIdx = inStream.findField("sample_id");
        int hammerColIdx = inStream.findField("hammer");
        int weightColIdx = inStream.findField("weight");
        // We need some counters.
        int hammersIn = 0;
        int badSamples = 0;
        // Loop through the input lines.
        log.info("Scanning {} for hammer data.", this.inFile);
        long lastMsg = System.currentTimeMillis();
        for (var line : inStream) {
            hammersIn++;
            // Get the set of expected repgens.
            String sampleId = line.get(sampColIdx);
            Set<String> expected = this.blackMap.get(sampleId);
            if (expected == null) {
                badSamples++;
                expected = NULL_SET;
            }
            // Get the weight.
            double weight = line.getDouble(weightColIdx);
            // Get the hammer itself and find its descriptor.
            String hammer = line.get(hammerColIdx);
            HammerData hammerData = retVal.computeIfAbsent(hammer, x -> new HammerData());
            hammerData.update(weight, expected);
            long now = System.currentTimeMillis();
            if (now - lastMsg >= 5000) {
                log.info("{} hammers read, {} unique, {} bad sample IDs.", hammersIn, retVal.size(), badSamples);
                lastMsg = now;
            }
        }
        log.info("{} total hammers read, {} unique, {} bad sample IDs.", hammersIn, retVal.size(), badSamples);
        return retVal;
    }

    /**
     * Write the hammer data to the output file.
     *
     * @param tempStream	print writer for temporary output file
     * @param hammerMap		map of hammer strings to descriptors
     */
    private void writeTempFile(PrintWriter tempStream, Map<String, HammerData> hammerMap) {
        int linesOut = 0;
        long lastMsg = System.currentTimeMillis();
        log.info("Writing reduced hammer data to {}.", this.tempFile);
        // Write the output header.
        tempStream.println("hammer\tweight\trepgens");
        // Loop through the hammer map, writing output.
        for (var hammerEntry : hammerMap.entrySet()) {
            String hammer = hammerEntry.getKey();
            HammerData desc = hammerEntry.getValue();
            tempStream.println(hammer + "\t" + desc.getWeight() + "\t" + StringUtils.join(desc.getExpected(), ','));
            linesOut++;
            long now = System.currentTimeMillis();
            if (now - lastMsg >= 5000) {
                log.info("{} hammers output.", linesOut);
                lastMsg = now;
            }
        }
        log.info("{} total hammers output.", linesOut);
    }

    @Override
    protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
        try {
            int hammersIn = 0;
            int genomeClose = 0;
            int hammerClose = 0;
            log.info("Processing reduced hammer set in {}.", this.tempFile);
            long lastMsg = System.currentTimeMillis();
            // Remember the kmer size.
            int kSize = hammerDb.getKmerSize();
            // Write the output header.
            writer.println("hammer\tweight\tgenome_close\thammer_close");
            // Loop through the temporary file.
            try (TabbedLineReader tempStream = new TabbedLineReader(this.tempFile)) {
                for (var line : tempStream) {
                    hammersIn++;
                    String hammer = line.get(0);
                    double weight = line.getDouble(1);
                    String[] expected = StringUtils.split(line.get(2), ',');
                    // Create the close-hammer set for this hammer.
                    Set<String> closeHammers = this.getCloseHammers(hammer);
                    // Search for the close hammers in the hammer database.
                    Iterator<String> iter = closeHammers.iterator();
                    boolean hFound = false;
                    while (! hFound && iter.hasNext())
                        hFound = (hammerDb.getSource(iter.next()) != null);
                    // Now "hFound" will be TRUE if this hammer is close to another hammer.  If it is, then it is axiomatically
                    // also close to one of the other repgens, since the hammer is IN another repgen.
                    if (hFound) {
                        hammerClose++;
                        writer.println(hammer + "\t" + weight + "\t1\t1");
                    } else {
                        // Find out if it is close to one of the expected genomes.
                        boolean gFound = false;
                        for (int i = 0; ! gFound && i < expected.length; i++) {
                            // Get and cache this expected genome.
                            Genome iGenome = this.repGtos.computeIfAbsent(expected[i], x -> this.repGenomes.getGenome(x));
                            // Loop through the contigs.
                            Iterator<Contig> contigIter = iGenome.getContigs().iterator();
                            while (! gFound && contigIter.hasNext()) {
                                Contig contig = contigIter.next();
                                // Try to find one of the hammers in the forward direction.
                                gFound = this.checkForKmers(contig.getSequence(), kSize, closeHammers);
                                if (! gFound) {
                                    // No luck, check the reverse direction.
                                    gFound = this.checkForKmers(contig.getRSequence(), kSize, closeHammers);
                                }
                            }
                        }
                        if (gFound) {
                            // Here we are genome-close but NOT hammer-close
                            writer.println(hammer + "\t" + weight + "\t1\t");
                            genomeClose++;
                        }
                    }
                    long now = System.currentTimeMillis();
                    if (now - lastMsg >= 5000) {
                        log.info("{} hammers processed, {} hammer-close, {} genome-close only.", hammersIn, hammerClose, genomeClose);
                        lastMsg = now;
                    }
                }
            }
        } finally {
            if (! this.keepFlag) {
                log.info("Deleting temporary file {}.", this.tempFile);
                FileUtils.forceDelete(this.tempFile);
            }
        }
    }

    /**
     * Compute the hammers close to the input hammer.
     *
     * @param hammer	hammer of interest
     *
     * @return a set of hammers that differ from the incoming hammer at exactly one position
     */
    private Set<String> getCloseHammers(String hammer) {
        final int kSize = hammer.length();
        // For each position, we get 1 close-hammer per additional nucleic acid.  This allocates
        // extra space for hashing.
        Set<String> retVal = new HashSet<String>(kSize * (4 * ACIDS.length) / 3 + 1);
        // We will create the close-hammers in this buffer.
        StringBuilder buffer = new StringBuilder(hammer);
        // Loop through the hammer, creating the close-hammers for each position.
        for (int i = 0; i < kSize; i++) {
            // Get the current character.
            char curr = hammer.charAt(i);
            for (char acid : ACIDS) {
                if (acid != curr) {
                    // Create a close-hammer.
                    buffer.setCharAt(i, acid);
                    retVal.add(buffer.toString());
                }
            }
            // Restore the buffer.
            buffer.setCharAt(i, curr);
        }
        return retVal;
    }

    /**
     * Determine whether or not any of the kmers in the specified set occur in the specified sequence.
     *
     * @param sequence		sequence to examine
     * @param kSize			kmer size to use
     * @param closeHammers	set of kmers to find
     *
     * @return TRUE if one of the kmers occurs at least once, else FALSE
     */
    private boolean checkForKmers(String sequence, int kSize, Set<String> closeHammers) {
        SequenceKmerIterable seqKmers = new SequenceKmerIterable(sequence, kSize);
        boolean retVal = seqKmers.stream().anyMatch(x -> closeHammers.contains(x));
        return retVal;
    }


}
