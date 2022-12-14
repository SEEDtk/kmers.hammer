/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command uses a protein-finder to compute hammers for a repgen set.  The standard input should contain the repgen set's
 * list file, which specifies the genomes in each neighborhood of a representative.  For each hammer role, this command
 * searches through that role's finder file for kmers that occur in most neighbors of a representative (worthiness) and few
 * genomes outside the neighborhood (precision).
 *
 * The hammers will be written to the standard output, in the standard two-column format for a hammer database.
 *
 * The positional parameter is the name of the protein-finder directory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing the repgen list data (if not STDIN)
 * -o	output file for the hammer list (if not STDOUT)
 * -K	kmer size for a hammer (default 20)
 *
 * --minWorth	minimum worthiness fraction for an acceptable hammer (default 0.50)
 * --minPrec	minimum precision fraction for an acceptable hammer (default 0.90)
 * --para		if specified, parallel processing will be used
 *
 * @author Bruce Parrello
 *
 */
public class HammerFinderProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerFinderProcessor.class);
    /** protein finder */
    private ProteinFinder finder;
    /** map of repgen genomes to neighbor sets */
    private Map<String, Set<String>> neighborMap;
    /** total number of genomes in the finder */
    private int genomeCount;
    /** ID of the current role */
    private String roleId;

    // COMMAND-LINE OPTIONS

    /** hammer size, in base pairs */
    @Option(name = "--kSize", aliases = { "-K", "--kmer" }, metaVar = "21", usage = "hammer DNA kmer size")
    private int kmerSize;

    /** minimum worthiness fraction for an acceptable hammer */
    @Option(name = "--minWorth", metaVar = "0.80", usage = "minimum worthiness fraction for an acceptable hammer")
    private double minWorth;

    /** minimum precision fraction for an acceptable hammer */
    @Option(name = "--minPrec", metaVar = "0.75", usage = "minimum precision fraction for an acceptable hammer")
    private double minPrec;

    /** TRUE to do parallel processing */
    @Option(name = "--para", usage = "if specified, parallel processing will be used")
    private boolean paraMode;

    /** protein-finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory")
    private File finderDir;

    /**
     * This object stores the scores for a hammer.
     */
    public static class Scores {
        /** worthiness score */
        private double worthiness;
        /** precision score */
        private double precision;

        /**
         * Construct a hammer score
         *
         * @param worth		worthiness ratio
         * @param prec		precision ratio
         */
        public Scores(double worth, double prec) {
            this.worthiness = worth;
            this.precision = prec;
        }

        /**
         * @return the worthiness
         */
        public double getWorthiness() {
            return this.worthiness;
        }

        /**
         * @return the precision
         */
        public double getPrecision() {
            return this.precision;
        }

    }
    @Override
    protected void setPipeDefaults() {
        this.kmerSize = 20;
        this.minPrec = 0.9;
        this.minWorth = 0.5;
        this.paraMode = false;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer size must be at least 2.");
        if (this.minPrec < 0.0 || this.minPrec > 1.0)
            throw new ParseFailureException("Minimum precision must be between 0 and 1.");
        if (this.minWorth < 0.0 || this.minWorth > 1.0)
            throw new ParseFailureException("Minimum worthiness must be between 0 and 1.");
        if (! this.finderDir.isDirectory())
            throw new FileNotFoundException("Finder directory " + this.finderDir + " is not found or invalid.");
        log.info("Loading protein-finder from {}.", this.finderDir);
        this.finder = new ProteinFinder(this.finderDir);
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Here we must read the representative list file.  For each representative genome, we create a set of its
        // neighbors.  To do this, we need to find two key input fields.
        int genomeColIdx = inputStream.findField("genome_id");
        int repColIdx = inputStream.findField("rep_id");
        log.info("Scanning repgen list file.");
        this.neighborMap = new HashMap<String, Set<String>>(1000);
        this.genomeCount = 0;
        // Connect each genome to its representative.
        for (var line : inputStream) {
            String genomeId = line.get(genomeColIdx);
            String repId = line.get(repColIdx);
            Set<String> neighborSet = this.neighborMap.computeIfAbsent(repId, x -> new HashSet<String>());
            neighborSet.add(genomeId);
            this.genomeCount++;
        }
        log.info("{} genomes sorted into {} neighborhoods.", this.genomeCount, this.neighborMap.size());
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Write the output header.
        writer.println("hammer\tfid\tworthiness\tprecision");
        // Now we are ready to begin.  For each role, we load all the sequences into memory, grouped
        // by representative genome.
        var fastaFileMap = this.finder.getFastas();
        for (var fastaEntry : fastaFileMap.entrySet()) {
            this.roleId = fastaEntry.getKey();
            File fastaFile = fastaEntry.getValue();
            log.info("Processing hammers for role {} using {}.", this.roleId, fastaFile);
            Map<String, List<HammerKmers>> genomeKmerMap = this.readFastaFile(fastaFile);
            // We process each representative genome separately, producing hammers from it.
            Stream<String> genomeStream;
            if (this.paraMode)
                genomeStream = this.neighborMap.keySet().parallelStream();
            else
                genomeStream = this.neighborMap.keySet().stream();
            genomeStream.forEach(x -> this.processGenomeHammers(x, genomeKmerMap, writer));
            log.info("Processing complete for role {}.", this.roleId);
        }
    }

    /**
     * Search for hammers that are worthy and precise in the specified genome.
     *
     * @param genomeId			ID of the representative genome for which we want kmers
     * @param genomeKmerMap		map of genome IDs to hammer-kmer objects
     * @param writer			output writer for hammers
     */
    private void processGenomeHammers(String genomeId, Map<String, List<HammerKmers>> genomeKmerMap, PrintWriter writer) {
        long start = System.currentTimeMillis();
        log.info("Searching for {} hammers in {}.", this.roleId, genomeId);
        // Get the kmers for the representative genome.
        var kmerList = genomeKmerMap.get(genomeId);
        if (kmerList != null) {
            // Get some counters.
            int unworthy = 0;
            int imprecise = 0;
            int processed = 0;
            // Get the set of neighbor genomes.
            Set<String> neighbors = this.neighborMap.get(genomeId);
            // Compute the denominator for the precision fraction.
            int foreignerCount = this.genomeCount - neighbors.size();
            // This will be the output set of kmers.
            Map<String, Scores> hammersFound = new HashMap<String, Scores>(10000);
            // This tracks the kmers we've already seen so we don't process them twice.
            Set<String> hammersSeen = new HashSet<String>(10000);
            // Loop through the feature-kmer list.
            for (HammerKmers kmers : kmerList) {
                for (String hammer : kmers.getKmers()) {
                    if (! hammersSeen.contains(hammer)) {
                        // Check this hammer for worthiness.
                        int worthCount = this.checkWorth(hammer, neighbors, genomeKmerMap);
                        double worthiness = worthCount / (double) neighbors.size();
                        if (worthiness < this.minWorth)
                            unworthy++;
                        else {
                            // The hammer is worthy.  Check for precision.
                            int errorCount = this.checkPrecision(hammer, neighbors, genomeKmerMap);
                            double precision = (foreignerCount - errorCount) / (double) foreignerCount;
                            if (precision < this.minPrec)
                                imprecise++;
                            else
                                hammersFound.put(hammer, new Scores(worthiness, precision));
                        }
                        // Insure we don't reprocess this hammer.
                        hammersSeen.add(hammer);
                        processed++;
                        if (log.isInfoEnabled() && processed % 300 == 0)
                            log.info("{} hammers processed for {} in {}.", processed, this.roleId, genomeId);
                    }
                }
                this.writeHammers(writer, kmers.getFid(), hammersFound);
            }
            log.info("{} {}-hammers found in {}. {} imprecise, {} unworthy.", hammersFound.size(), this.roleId,
                    genomeId, imprecise, unworthy);
        }
        if (log.isInfoEnabled()) {
            Duration d = Duration.ofMillis(System.currentTimeMillis() - start);
            log.info("{} to process {} for {}.", d, genomeId, this.roleId);
        }

    }

    /**
     * Count the number of times the specified hammer is found in a neighbor genome.
     *
     * @param hammer		hammer to look for
     * @param neighbors		set of neighbor genome IDs
     * @param genomeKmerMap	map of genome IDs to kmer objects
     *
     * @return the number of neighbors that contained the hammer
     */
    private int checkWorth(String hammer, Set<String> neighbors, Map<String, List<HammerKmers>> genomeKmerMap) {
        // This is a fairly complex stream.  For each neighbor, we get its list of hammer kmers (which is almost
        // always a singleton), using an empty list if the neighbor does not have a sequence.  If we find the
        // hammer in a genome's kmers, we count it.  The worthiness count is the number of neighbor genomes
        // containing the hammer.  The representative genome is included, so the count is always at least 1.
        int retVal = (int) neighbors.stream().map(x -> genomeKmerMap.getOrDefault(x, Collections.emptyList()))
                .filter(x -> this.checkHammer(hammer, x, HammerKmers.Mode.WORTHINESS)).count();
        // Return the count.
        return retVal;
    }

    /**
     * Count the number of times the specified hammer is found in a non-neighbor genome.
     *
     * @param hammer		hammer to look for
     * @param neighbors		set of neighbor genomes to skip
     * @param genomeKmerMap	map of genome IDs to kmer objects
     *
     * @return the number of non-neighbors that contained the hammer
     */
    private int checkPrecision(String hammer, Set<String> neighbors, Map<String, List<HammerKmers>> genomeKmerMap) {
        int retVal = 0;
        // Loop through the neighborhoods other than the representative's neighborhood.
        for (var genomeKmerEntry : genomeKmerMap.entrySet()) {
            // Only proceed if this is a non-neighbor.
            if (! neighbors.contains(genomeKmerEntry.getKey())) {
                // Check the kmers for this non-neighbor genome.
                if (this.checkHammer(hammer, genomeKmerEntry.getValue(), HammerKmers.Mode.PRECISION))
                    retVal++;
            }
        }
        return retVal;
    }

    /**
     * Check a list of hammer-kmer objects for a particular hammer.
     *
     * @param hammer		hammer to check for
     * @param kmerList		list of hammer-kmer objects to check
     * @param type			type of check (WORTHINESS or PRECISION)
     *
     * @return TRUE if the hammer was found, else FALSE
     */
    private boolean checkHammer(String hammer, Collection<HammerKmers> kmerList, HammerKmers.Mode type) {
        boolean retVal = kmerList.stream().anyMatch(x -> type.check(hammer, x));
        return retVal;
    }

    /**
     * Write all the hammers found to the output.  This method is single-threaded because the output is
     * not a thread-safe resource.
     *
     * @param writer		output writer
     * @param fid			ID of the feature containing the hammer
     * @param hammersFound	set of hammers found with their scores
     */
    private synchronized void writeHammers(PrintWriter writer, String fid, Map<String, Scores> hammersFound) {
        for (var hammerInfo : hammersFound.entrySet()) {
            String hammer = hammerInfo.getKey();
            Scores scoreInfo = hammerInfo.getValue();
            writer.println(hammer + "\t" + fid + "\t" + Double.toString(scoreInfo.getWorthiness())
                    + "\t" + Double.toString(scoreInfo.getPrecision()));
        }
        writer.flush();
    }

    /**
     * Read all the sequences from a FASTA file into memory and organize them by genome ID.
     *
     * @param fastaFile		FASTA file containing the DNA sequences for a role
     *
     * @return a map from each genome ID to the list of kmer sets of its sequences
     *
     * @throws IOException
     */
    private Map<String, List<HammerKmers>> readFastaFile(File fastaFile) throws IOException {
        var retVal = new HashMap<String, List<HammerKmers>>(this.genomeCount * 4 / 3 + 1);
        log.info("Reading role hammers from {}.", fastaFile);
        int hCount = 0;
        try (FastaInputStream inStream = new FastaInputStream(fastaFile)) {
            for (Sequence seq : inStream) {
                String fid = seq.getLabel();
                HammerKmers kmerSet = new HammerKmers(fid, seq.getSequence());
                String genomeId = Feature.genomeOf(fid);
                // Add the hammer-kmer object to the genome's sequence list.  Most of the time there will only be
                // one, but we allow for more, hence the use of a linked list.
                retVal.computeIfAbsent(genomeId, x -> new LinkedList<HammerKmers>()).add(kmerSet);
                hCount++;
            }
        }
        log.info("{} sequences for {} genomes read from {}.", hCount, retVal.size(), fastaFile);
        return retVal;
    }

}
