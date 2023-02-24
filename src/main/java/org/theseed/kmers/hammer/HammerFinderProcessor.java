/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
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
 * --temp		name of a temporary directory for working files (default "Temp" in the current directory)
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
    /** hammer scoring map */
    private Map<String, HammerScore> hammerMap;
    /** two-level map from genomeId -> featureID -> feature DNA for all sequences of the current role */
    private Map<String, Map<String, String>> sequenceMap;
    /** total number of genomes in the finder */
    private int genomeCount;
    /** ID of the current role */
    private String roleId;

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

    /** protein-finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory")
    private File finderDir;


    @Override
    protected void setPipeDefaults() {
        this.kmerSize = 20;
        this.minPrec = 0.9;
        this.minWorth = 0.0;
        this.minPeers = 1;
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
        // Insure the hammer scoring system has the total genome count.
        HammerScore.setTotalGenomes(this.genomeCount);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Write the output header.
        writer.println("hammer\tfid\tstrength\tprecision\thits");
        // Now we are ready to begin.  For each role, we load all the sequences into memory, grouped
        // by representative genome.
        var fastaFileMap = this.finder.getFastas();
        for (var fastaEntry : fastaFileMap.entrySet()) {
            this.roleId = fastaEntry.getKey();
            File fastaFile = fastaEntry.getValue();
            log.info("Processing hammers for role {} using {}.", this.roleId, fastaFile);
            // Read all the sequences from the FASTA file into the sequence map.
            this.sequenceMap = this.readSequences(fastaFile);
            // Generate the potential hammers.
            this.hammerMap = this.getPotentialHammers();
            // Compute worthiness and precision.
            this.analyzeHammers();
            // Write the worthwhile hammers.
            this.writeResults(writer);
        }
    }

    /**
     * Read the specified finder file and return all the sequences organized by genome.
     *
     * @param fastaFile		FASTA file of DNA sequences for the current role
     *
     * @return a map of feature IDs to sequences for each genome
     *
     * @throws IOException
     */
    private Map<String, Map<String, String>> readSequences(File fastaFile) throws IOException {
        // Create the return map.  There will be one entry per genome in the system.
        var retVal = new HashMap<String, Map<String, String>>(this.genomeCount * 4 / 3 + 1);
        // Loop through the FASTA file.
        log.info("Reading sequences from {} for role {}.", fastaFile, this.roleId);
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
                if (log.isInfoEnabled() && inCount % 5000 == 0)
                    log.info("{} sequences processed for {} genomes.", inCount, retVal.size());
            }
        }
        return retVal;
    }

    /**
     * Compute the potential hammers.  We do this by scanning all the sequences belonging to one of the representative
     * genomes (that is, the ones whose IDs are in the key set of the neighborhood map).  If a sequence is found in
     * more than one such genome, it is removed from the map.  This is the method that provides the greatest memory
     * strain.  All subsequent methods merely reduce the size of the hammer map.
     *
     * @return a map from potential hammers to scoring objects
     */
    private Map<String, HammerScore> getPotentialHammers() {
        // The potential hammers are stored in here.
        Map<String, HammerScore> retVal = new HashMap<String, HammerScore>(this.neighborMap.size() * 4000);
        // Common kmers are kept in here.
        Set<String> commonSet = new HashSet<String>(this.neighborMap.size() * 150);
        // Loop through the representatives.
        for (var repEntry : this.neighborMap.entrySet()) {
            String repId = repEntry.getKey();
            var repSequenceMap = this.sequenceMap.get(repId);
            if (repSequenceMap != null) {
                // Get the neighbor set.  We store it in the score objects.
                Set<String> neighbors = repEntry.getValue();
                log.info("Scanning representative genome {} with {} neighbors.", repId, neighbors.size() - 1);
                // Loop through all the sequences for this representative.  We expect one, maybe two.
                // Note that we don't check for ambiguity characters, since the finder is built without them.
                int commonCount = 0;
                int conflictCount = 0;
                int kmerCount = 0;
                for (var seqEntry : repSequenceMap.entrySet()) {
                    String fid = seqEntry.getKey();
                    KmerSeries kmers = new KmerSeries(seqEntry.getValue(), this.kmerSize);
                    for (String kmer : kmers) {
                        kmerCount++;
                        if (commonSet.contains(kmer)) {
                            // Here the kmer has already been identified as bad.
                            commonCount++;
                        } else {
                            HammerScore score = retVal.get(kmer);
                            if (score == null) {
                                // Here the kmer is new.  Add it to the return map.
                                retVal.put(kmer, new HammerScore(fid, neighbors));
                            } else if (score.isDisqualifyingHit(fid)) {
                                // Here the kmer has been seen before in a different genome.
                                retVal.remove(kmer);
                                commonSet.add(kmer);
                                conflictCount++;
                            }
                        }
                    }
                }
                log.info("{} kmers processed:  {} were common, {} conflicts found.", kmerCount, commonCount, conflictCount);
                // Now delete the sequences for this rep genome.  We don't need them again.
                this.sequenceMap.remove(repId);
            }
        }
        log.info("{} potential kmers found for {}.", retVal.size(), this.roleId);
        return retVal;
    }

    /**
     * Now we score the hammers.  We run through ALL the sequences in genome order, updating the scores.
     */
    private void analyzeHammers() {
        log.info("Scoring potential hammers for role {}.", this.roleId);
        int seqCount = 0;
        for (var genomeEntry : this.sequenceMap.entrySet()) {
            String genomeId = genomeEntry.getKey();
            for (String seq : genomeEntry.getValue().values()) {
                KmerSeries kmers = new KmerSeries(seq, this.kmerSize);
                for (String kmer : kmers) {
                    HammerScore score = this.hammerMap.get(kmer);
                    if (score != null)
                        score.recordHit(genomeId);
                }
                seqCount++;
                if (log.isInfoEnabled() && seqCount % 5000 == 0)
                    log.info("{} sequences processed for role {}.", seqCount, this.roleId);
            }
        }
    }

    /**
     * Here we output the hammers with their scores.
     *
     * @param writer	output writer for the report
     */
    private void writeResults(PrintWriter writer) {
        int worthless = 0;
        int imprecise = 0;
        int kept = 0;
        int lowSize = 0;
        for (var hammerEntry : this.hammerMap.entrySet()) {
            String hammer = hammerEntry.getKey();
            HammerScore score = hammerEntry.getValue();
            String fid = score.getFid();
            String genomeId = Feature.genomeOf(fid);
            // We must insure the representative has a good neighborhood.
            if (this.neighborMap.get(genomeId).size() < this.minPeers)
                lowSize++;
            else {
                double worth = score.getWorthiness();
                double strength = score.getStrength();
                double prec = score.getPrecision();
                if (worth < this.minWorth)
                    worthless++;
                else if (prec < this.minPrec)
                    imprecise++;
                else {
                    kept++;
                    writer.println(hammer + "\t" + score.getFid() + "\t" + Double.toString(strength)
                            + "\t" + Double.toString(prec) + "\t" + Integer.toString(score.getGoodHits()));
                    if (log.isInfoEnabled() && kept % 30000 == 0)
                        log.info("{} hammers output, {} unworthy, {} imprecise.", kept, worthless, imprecise);
                }
            }
        }
        log.info("{} hammers output, {} unworthy, {} imprecise, {} rejected for insufficient neighborhood.",
                kept, worthless, imprecise, lowSize);
    }

}
