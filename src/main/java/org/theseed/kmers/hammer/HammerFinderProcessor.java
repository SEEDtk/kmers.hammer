/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerMap;
import org.theseed.proteins.hammer.HammerScore;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.KmerSeries;
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
 * The hammers will be written to the standard output, in the standard three-column format for a hammer database.
 *
 * The hammers are generated one role at a time.  For each role, we read in all the finder sequences and sort them by
 * genome.  We then build a map from hammers in the representative genomes to scoring objects.
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
 * --minWorth	minimum worthiness fraction for an acceptable hammer (default 0.00)
 * --minPrec	minimum precision fraction for an acceptable hammer (default 0.90)
 * --minSize	minimum size for a neighborhood to be hammer-worthy (default 1)
 * --maxRepeat	maximum percent of a hammer that can belong to a single base pair type (default 0.70)
 * --para		if specified, parallel processing will be used, which increases memory size but improves performance
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
    /** maximum repeat count */
    private int maxRepeatCount;
    /** number of hammers scanned from representative genomes */
    private long repScanCount;
    /** number of hammers rejected due to low complexity */
    private long rejectCount;
    /** number of hammers common between representatives */
    private long conflictCount;
    /** number of kmers scanned in neighbor genomes */
    private long kmerCount;
    /** hammer map used to accumulate potential hammers */
    private HammerMap<HammerScore> hammerMap;

    // COMMAND-LINE OPTIONS

    /** hammer size, in base pairs */
    @Option(name = "--kSize", aliases = { "-K", "--kmer" }, metaVar = "21", usage = "hammer DNA kmer size")
    private int kmerSize;

    /** minimum worthiness fraction for an acceptable hammer */
    @Option(name = "--minWorth", metaVar = "0.95", usage = "minimum worthiness fraction for an acceptable hammer")
    private double minWorth;

    /** minimum precision fraction for an acceptable hammer */
    @Option(name = "--minPrec", metaVar = "1.0", usage = "minimum precision fraction for an acceptable hammer")
    private double minPrec;

    /** minimum size of a neighborhood for hammers to be worthwhile */
    @Option(name = "--minSize", metaVar = "50", usage = "minimum neighborhood size for a representative")
    private int minPeers;

    /** maximum fraction of a hammer allowed for a particular nucleotide */
    @Option(name = "--maxRepeat", metaVar = "0.60", usage = "maximum fraction of a hammer allowed for any one nucleotide type")
    private double maxRepeat;

    /** if specified, parallel processing will be used */
    @Option(name = "--para", usage = "if specified, parallel processing will be used")
    private boolean paraFlag;

    /** protein-finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory")
    private File finderDir;


    @Override
    protected void setPipeDefaults() {
        this.kmerSize = 20;
        this.minPrec = 0.9;
        this.minWorth = 0.0;
        this.minPeers = 1;
        this.maxRepeat = 0.70;
        this.paraFlag = false;
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
        if (this.minPeers <= 0)
            throw new ParseFailureException("Minimum neighborhood size must be positive.");
        if (this.maxRepeat <= 0.25 || this.maxRepeat > 1.0)
            throw new ParseFailureException("Maximum repeat fraction must be between 0.25 and 1.0.");
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
        HammerScore.setTotalGenomes(this.genomeCount);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Compute the repeat-count limit.  A value greater than this will cause the hammer to be
        // discarded.
        this.maxRepeatCount = (int) (this.maxRepeat * this.kmerSize);
        // Initialize the counters.
        this.repScanCount = 0;
        this.rejectCount = 0;
        this.conflictCount = 0;
        this.kmerCount = 0;
        // Create the master hammer map.
        this.hammerMap = new HammerMap<HammerScore>(this.kmerSize);
        // Now we are ready to begin.  For each role, we load all the sequences into memory, grouped
        // by representative genome.
        var fastaFileMap = this.finder.getFastas();
        var fastaStream = fastaFileMap.entrySet().parallelStream();
        if (this.paraFlag)
            fastaStream = fastaStream.parallel();
        fastaStream.forEach(x -> this.processRole(x));
        log.info("{} hammers scanned in {} representative genomes.  {} rejected due to low complexity, {} due to conflicts, {} kept.",
                this.repScanCount, this.neighborMap.size(), this.rejectCount, this.conflictCount, this.hammerMap.size());
        log.info("{} total kmers scanned in neighbor genomes.", this.kmerCount);
        log.info("Hammer map has a load factor of {} and overload factor of {}.", this.hammerMap.loadFactor(),
                this.hammerMap.overloadFactor());
        // Now we write the results.
        this.writeHammers(writer);
    }

    /**
     * Process the hammers for a single role.
     *
     * @param fastaEntry	map entry containing the role ID and the role's finder FASTA file
     *
     */
    private void processRole(Entry<String, File> fastaEntry) {
        String roleId = fastaEntry.getKey();
        File fastaFile = fastaEntry.getValue();
        log.info("Processing hammers for role {} using {}.", roleId, fastaFile);
        // Read all the sequences from the FASTA file into the sequence map.  We need to convert the IO
        // exception so we can use this method in a stream.
        Map<String, Map<String, String>> sequenceMap;
        try {
            sequenceMap = this.readSequences(fastaFile, roleId);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        // Now we have all the DNA sequences for this role in the sequence map, organized by genome ID.
        // We run through all the repgens, storing the hammers in the map.  Then we run through all
        // the others, counting the hits.  First, the repgens.
        log.info("Scanning hammers for role {} in representative genomes.", roleId);
        this.findHammers(roleId, sequenceMap);
        // Now, the other genomes.
        log.info("Counting hammer occurrences for role {} in non-representative genomes.", roleId);
        this.countHammers(roleId, sequenceMap);
    }

    /**
     * Read the specified finder file and return all the sequences organized by genome.
     *
     * @param fastaFile		FASTA file of DNA sequences for the current role
     * @param roleId		ID of the current role
     *
     * @return a map of feature IDs to sequences for each genome
     *
     * @throws IOException
     */
    private Map<String, Map<String, String>> readSequences(File fastaFile, String roleId) throws IOException {
        long lastMsg = System.currentTimeMillis();
        // Create the return map.  There will be one entry per genome in the system.
        var retVal = new HashMap<String, Map<String, String>>(this.genomeCount * 4 / 3 + 1);
        // Loop through the FASTA file.
        log.info("Reading sequences from {} for role {}.", fastaFile, roleId);
        try (var inStream = new FastaInputStream(fastaFile)) {
            int inCount = 0;
            for (Sequence seq : inStream) {
                // Get the genome ID for this sequence.
                String fid = seq.getLabel();
                String genomeId = Feature.genomeOf(fid);
                // Get the genome's sequence map.
                Map<String,String> gMap = retVal.computeIfAbsent(genomeId, x -> new TreeMap<String, String>());
                gMap.put(fid, seq.getSequence());
                inCount++;
                if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 5000) {
                    log.info("{} sequences for role {} read for {} genomes.", inCount, roleId, retVal.size());
                    lastMsg = System.currentTimeMillis();
                }
            }
        }
        return retVal;
    }

    /**
     * Find all the hammers in the representative genomes for the current finder.  A hammer that is already in the map for
     * a different genome will automatically be marked bad.  Sequences that fail the complexity test will be skipped.
     * There is no need to scan for ambiguity characters because the finder is built without them.
     *
     * @param roleId		ID of the role being processed
     * @param sequenceMap	two-level map from genome ID -> feature ID -> DNA sequence
     */
    private void findHammers(String roleId, Map<String, Map<String, String>> sequenceMap) {
        // Create some counters.
        long scanCount = 0;
        long badCount = 0;
        long commonCount = 0;
        // Loop through the repgen genomes.  We need the genome ID and the neighbor set for each.
        for (var repEntry : this.neighborMap.entrySet()) {
            String repId = repEntry.getKey();
            Set<String> neighbors = repEntry.getValue();
            // If the neighborhood is too small, we have to remember all hammers as bad.
            boolean alwaysBad = (neighbors.size() < this.minPeers);
            // Does this repgen have sequences for this protein?
            var fidMap = sequenceMap.get(repId);
            if (fidMap != null) {
                // Yes, it does. In almost every case, there will be a single feature in the finder for this genome.
                // If there is more than one, we default to the first.
                var fid = fidMap.keySet().iterator().next();
                // Get all the kmers for this protein's sequences.
                Set<String> kmers = this.getHammers(fidMap.values());
                // For each kmer, if it is new, we add it as a hammer.  If it is old, we mark it as invalid.
                // We must, however, check for low-complexity sequences.
                for (String kmer : kmers) {
                    scanCount++;
                    int bestBaseCount = HammerMap.commonBaseCount(kmer);
                    if (bestBaseCount > this.maxRepeatCount)
                        badCount++;
                    else {
                        // Here the hamemr is valid.  Check it against the map.  If it already exists, remember it
                        // as bad.  Otherwise, add it to the map.  Map updates are thread-safe, and even if we bad-flag
                        // the score, it is the only operation allowed on that flag, so thread safety does not matter.
                        boolean newHammer = this.hammerMap.update(kmer, x -> x.setBadHammer(),
                                x -> new HammerScore(fid, roleId, neighbors, alwaysBad));
                        if (!newHammer)
                            commonCount++;
                    }
                }
            }
        }
        log.info("{} kmers scanned, {} rejected, {} in common for role {}.", scanCount, badCount, commonCount, roleId);
        // Update the counters.
        synchronized (this) {
            this.repScanCount += scanCount;
            this.rejectCount += badCount;
            this.conflictCount += commonCount;
        }
    }

    /**
     * Count all the hammers in the represented genomes so we can determine the scores.
     *
     * @param roleId	ID of the role being processed
     * @param sequenceMap	two-level map from genome ID -> feature ID -> DNA sequence
     */
    private void countHammers(String roleId, Map<String, Map<String, String>> sequenceMap) {
        // This timer is used to insure we get a message every 5 seconds.
        long lastMsg = System.currentTimeMillis();
        // Loop through all the genomes, counting the hammers in the non-representatives.
        long kCount = 0;
        for (var genomeEntry : sequenceMap.entrySet()) {
            String genomeId = genomeEntry.getKey();
            if (! this.neighborMap.containsKey(genomeId)) {
                // Here we are definitely not a representative genome.
                Set<String> kmers = this.getHammers(genomeEntry.getValue().values());
                for (String kmer : kmers) {
                    kCount++;
                    HammerScore score = this.hammerMap.get(kmer);
                    if (score != null)
                        score.recordHit(genomeId, roleId);
                }
                if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 5000) {
                    log.info("{} kmers scanned for role {}.", kCount, roleId);
                    lastMsg = System.currentTimeMillis();
                }
            }
        }
        log.info("Total for role {}: {} kmers scanned from non-representative genomes.", roleId, kCount);
        // Update the counter.
        synchronized (this) {
            this.kmerCount += kCount;
        }
    }

    /**
     * @return the set of hammers in a group of sequences
     *
     * @param seq	group of sequences to scan
     */
    private Set<String> getHammers(Collection<String> seqs) {
        Set<String> retVal = new HashSet<String>(seqs.size() * 1400);
        KmerSeries kmers = new KmerSeries(seqs, this.kmerSize);
        for (String kmer : kmers)
            retVal.add(kmer);
        return retVal;
    }

    /**
     * Write the hammers to the output.  This method happens only once, so unlike the rest, it is not thread-safe.
     *
     * @param writer	output file for the hammer list
     */
    private void writeHammers(PrintWriter writer) {
        log.info("Writing output (this will take many minutes).");
        long hammersIn = 0;
        long lowPrec = 0;
        long lowWorth = 0;
        long hammersOut = 0;
        // Write the output header.
        writer.println("hammer\tfid\tstrength\tprecision\tworth\thits\trole");
        // Loop through the hammers in the hammer map.
        for (var hammerEntry : this.hammerMap) {
            String hammer = hammerEntry.getKey();
            HammerScore score = hammerEntry.getValue();
            // Only proceed if the hammer is good.
            if (! score.isBadHammer()) {
                double worth = score.getWorthiness();
                double precision = score.getPrecision();
                String fid = score.getFid();
                hammersIn++;
                // Filter out the hammers that aren't useful.
                if (worth < this.minWorth)
                    lowWorth++;
                else if (precision < this.minPrec)
                    lowPrec++;
                else {
                    // Here we have a good hammer.
                    writer.println(hammer + "\t" + fid + "\t" + score.getStrength() + "\t"
                            + precision + "\t" + worth + "\t" + score.getGoodHits() + "\t"
                            + score.getRoleId());
                    hammersOut++;
                }
            }
        }
        log.info("{} hammers checked, {} rejected due to low precision, {} due to low worth, {} output.",
                hammersIn, lowPrec, lowWorth, hammersOut);
    }


}
