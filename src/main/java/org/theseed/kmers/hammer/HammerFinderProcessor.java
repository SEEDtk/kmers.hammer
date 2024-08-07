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
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ForkJoinPool;
import java.util.Set;
import java.util.TreeMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.hammer.HammerMap;
import org.theseed.proteins.hammer.HammerScore;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.KmerSeries;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.sequence.seeds.filters.HammerDupFilter;
import org.theseed.sequence.seeds.filters.HammerFeatureFilter;
import org.theseed.utils.BasePipeProcessor;

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
 * -t	type of genome source for the optional cleaning step (default DIR)
 *
 * --minWorth	minimum worthiness fraction for an acceptable hammer (default 0.00)
 * --minPrec	minimum precision fraction for an acceptable hammer (default 0.90)
 * --minSize	minimum size for a neighborhood to be hammer-worthy (default 1)
 * --maxRepeat	maximum percent of a hammer that can belong to a single base pair type (default 0.70)
 * --para		maximum number of parallel threads to use during cleaning
 * --clean		if specified, the name of a genome source for the representative genomes, to be used to clean the hammer set
 * --sType		type of strength computation to use (default RATIO_BASED)
 * --reject		if specified, a file to contain the hammers rejected due to conflicts
 * --pegFilter	filter for features to use in the finder (default ALL)
 * --dupFilter	algorithm to use for processing multi-occurring roles in a genome (default KEEP)
 * --maxCount	maximum number of times a hammer can be found in its own genome during cleaning (default 10)
 * --anchor		remove hammers that are very close (1 nucleotide difference) to other hammers
 *
 * @author Bruce Parrello
 *
 */
public class HammerFinderProcessor extends BasePipeProcessor implements HammerFeatureFilter.IParms, HammerDupFilter.IParms {

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
    /** genome source for representative genomes, or NULL if there is no cleaning step */
    private GenomeSource repGenomes;
    /** feature filter to use */
    private HammerFeatureFilter pegFilter;
    /** duplicate-role filter to use */
    private HammerDupFilter dupFilter;
    /** custom thread pool for parallel processing */
    private ForkJoinPool threadPool;

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
    @Option(name = "--para", usage = "number of threads to use during cleaning phase")
    private int maxThreads;

    /** if specified, the file or directory containing a representative genome source for cleaning */
    @Option(name = "--clean", metaVar = "GTOXXX", usage = "if specified the file or directory containing the representative genomes for cleaning")
    private File repDir;

    /** type of genome source for the optional cleaning directory */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source for optional cleaning genomes")
    private GenomeSource.Type sourceType;

    /** type of strength computation to use */
    @Option(name = "--sType", usage = "type of strength computation to use")
    private HammerScore.Type strengthType;

    /** optional rejected-hammer file */
    @Option(name = "--reject", usage = "optional output file for rejected hammers")
    private File rejectFile;

    /** hammer-feature filter used to select features for hammer analysis */
    @Option(name = "--pegFilter", usage = "filter for finder features to use for hammer generation")
    private HammerFeatureFilter.Type pegFilterType;

    /** hammer-feature filter used to handle duplicate features for a role */
    @Option(name = "--dupFilter", usage = "filter for duplicate features to select for hammer generation")
    private HammerDupFilter.Type dupFilterType;

    /** maximum number of times a hammer can duplicate in a genome */
    @Option(name = "--maxCount", metaVar = "1", usage = "maximum number of duplicates allowed for a hammer in a genome")
    private int maxCount;

    /** if specified, hammers that are very close to other hammers will be removed */
    @Option(name = "--anchor", usage = "if specified, hammers that are very close to other hammers are removed")
    private boolean anchorFlag;

    /** protein-finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory", required = true)
    private File finderDir;


    @Override
    protected void setPipeDefaults() {
        this.kmerSize = 20;
        this.minPrec = 0.9;
        this.minWorth = 0.0;
        this.minPeers = 1;
        this.maxRepeat = 0.70;
        this.maxThreads = Runtime.getRuntime().availableProcessors();
        this.repDir = null;
        this.sourceType = GenomeSource.Type.DIR;
        this.strengthType = HammerScore.Type.RATIO_BASED;
        this.rejectFile = null;
        this.pegFilterType = HammerFeatureFilter.Type.ALL;
        this.dupFilterType = HammerDupFilter.Type.KEEP;
        this.maxCount = 10;
        this.anchorFlag = false;
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
        if (this.maxCount < 1)
            throw new ParseFailureException("Maximum-duplicate count must be at least 1.");
        // Set up the cleaning directory.
        if (this.repDir == null) {
            log.info("No cleaning step will be performed.");
            this.repGenomes = null;
        } else {
            log.info("Hammers will be cleaned using representative genomes in {}.", this.repDir);
            this.repGenomes = this.sourceType.create(this.repDir);
        }
        // Initialize the hammer-feature filters.
        this.pegFilter = this.pegFilterType.create(this);
        this.dupFilter = this.dupFilterType.create(this);
        // Finally, set up the protein finder.
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
        // Create the custom thread pool.
        if (this.maxThreads == 1)
            this.threadPool = null;
        else {
            this.threadPool = new ForkJoinPool(this.maxThreads);
            log.info("Parallel processing selected with {} threads.", this.maxThreads);
        }
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
        // by representative genome, and put the good ones in the hammer map.
        var fastaFileMap = this.finder.getFastas();
        var fastaStream = fastaFileMap.entrySet().stream();
        if (this.threadPool != null)
            fastaStream = fastaStream.parallel();
        fastaStream.forEach(x -> this.findRoleHammers(x));
        log.info("{} hammers scanned in {} representative genomes.  {} rejected due to low complexity, {} due to conflicts, {} kept.",
                this.repScanCount, this.neighborMap.size(), this.rejectCount, this.conflictCount, this.hammerMap.size());
        // Now we run through it all again, scoring the hammers found.  We have to do this in separate runs to avoid a race
        // condition between finding and counting relating to bad hammers.
        fastaStream = fastaFileMap.entrySet().stream();
        if (this.threadPool != null)
            fastaStream = fastaStream.parallel();
        fastaStream.forEach(x -> this.countRoleHammers(x));
        log.info("{} total kmers scanned in neighbor genomes. Conflict count is now {}.", this.kmerCount, this.conflictCount);
        log.info("Hammer map has a load factor of {} and overload factor of {}.", this.hammerMap.loadFactor(),
                this.hammerMap.overloadFactor());
        if (this.anchorFlag) {
            log.info("Scanning to remove non-anchor hammers.");
            int removeCount = this.hammerMap.anchorize();
            log.info("{} non-anchor hammers marked.", removeCount);
        }
        if (this.repGenomes != null) {
            // Here we have to clean the hammers.  We scan each representative genome, and delete any hammer found outside
            // its target genome.
            var repIdStream = this.neighborMap.keySet().stream();
            int cleaned;
            if (this.threadPool == null)
                cleaned = repIdStream.mapToInt(x -> this.cleanHammers(x)).sum();
            else try {
                // Yes, use a parallel stream.
                cleaned = this.threadPool.submit(() ->
                        repIdStream.parallel().mapToInt(x -> this.cleanHammers(x)).sum()).get();
            } finally {
                this.threadPool.shutdown();
            }
            log.info("{} hammers were marked bad during cleaning.", cleaned);
            log.info("Final hammer map has a load factor of {} and overload factor of {}.", this.hammerMap.loadFactor(),
                    this.hammerMap.overloadFactor());
        }
        // Now we write the results.
        this.writeHammers(writer);
    }

    /**
     * Find the hammers for a single role.
     *
     * @param fastaEntry	map entry containing the role ID and the role's finder FASTA file
     *
     */
    private void findRoleHammers(Entry<String, File> fastaEntry) {
        String roleId = fastaEntry.getKey();
        File fastaFile = fastaEntry.getValue();
        // Read all the sequences from the FASTA file into the sequence map.
        Map<String, Map<String, String>> sequenceMap = this.readSequences(fastaFile, roleId);
        // Now we have all the DNA sequences for this role in the sequence map, organized by genome ID.
        // We run through all the repgens, storing the hammers in the map.  Then we run through all
        // the others, counting the hits.  First, the repgens.
        log.info("Scanning hammers for role {} in representative genomes.", roleId);
        this.findHammers(roleId, sequenceMap);
    }

    /**
     * Score the hammers for a single role.
     *
     * @param fastaEntry	map entry containing the role ID and the role's finder FASTA file
     *
     */
    private void countRoleHammers(Entry<String, File> fastaEntry) {
        String roleId = fastaEntry.getKey();
        File fastaFile = fastaEntry.getValue();
        // Read all the sequences from the FASTA file into the sequence map.
        Map<String, Map<String, String>> sequenceMap = this.readSequences(fastaFile, roleId);
        // Count the hammers in the non-representative genomes to score them.
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
     */
    private Map<String, Map<String, String>> readSequences(File fastaFile, String roleId) {
        // Create the return map.  There will be one entry per genome in the system.
        var retVal = new HashMap<String, Map<String, String>>(this.genomeCount * 4 / 3 + 1);
        // Loop through the FASTA file.
        log.info("Reading sequences from {} for role {}.", fastaFile, roleId);
        try (var inStream = new FastaInputStream(fastaFile)) {
            int inCount = 0;
            int filterCount = 0;
            for (Sequence seq : inStream) {
                // Apply the feature filter.
                boolean keep = this.pegFilter.check(seq);
                if (! keep)
                    filterCount++;
                else {
                    // Get the genome ID for this sequence.
                    String fid = seq.getLabel();
                    String genomeId = Feature.genomeOf(fid);
                    // Get the genome's sequence map.
                    Map<String,String> gMap = retVal.computeIfAbsent(genomeId, x -> new TreeMap<String, String>());
                    gMap.put(fid, seq.getSequence());
                }
                inCount++;
            }
            log.info("{} sequences read for role {}. {} removed by filter.", inCount, roleId, filterCount);
        } catch (IOException e) {
            // Convert the exception to unchecked so we can use this in a stream.
            throw new UncheckedIOException(e);
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
        long skipCount = 0;
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
                // If there is more than one, we only generate hammers for the ones chosen by the dup-filter.
                Map<String, String> filtered = this.dupFilter.filter(fidMap);
                // Count the filtered features.
                int newCount = filtered.size();
                skipCount += fidMap.size() - newCount;
                // Only proceed if any are left.
                if (newCount > 0) {
                    var fid = fidMap.keySet().iterator().next();
                    // Get all the kmers for this protein's sequences.
                    Set<String> kmers = this.getHammers(filtered.values());
                    // For each kmer, if it is new, we add it as a hammer.  If it is old, we mark it as invalid.
                    // We must, however, check for low-complexity sequences.
                    for (String kmer : kmers) {
                        scanCount++;
                        int bestBaseCount = HammerMap.commonBaseCount(kmer);
                        if (bestBaseCount > this.maxRepeatCount)
                            badCount++;
                        else {
                            // Here the hammer is valid.  Check it against the map.  If it already exists, remember it
                            // as bad.  Otherwise, add it to the map.  The update operation is thread-safe.
                            boolean newHammer = this.hammerMap.update(kmer, x -> x.setBadHammer(),
                                    x -> this.strengthType.create(fid, roleId, neighbors, alwaysBad));
                            if (!newHammer)
                                commonCount++;
                        }
                    }
                }
            }
        }
        log.info("{} kmers scanned, {} rejected, {} in common for role {}, {} duplicate features skipped.", scanCount, badCount, commonCount, roleId, skipCount);
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
                if (log.isInfoEnabled() && System.currentTimeMillis() - lastMsg >= 10000) {
                    log.info("{} kmers counted for role {}.", kCount, roleId);
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
     * This method scans a genome for hammers.  A hammer found in a genome that is different from its target
     * genome is removed.  A hammer that occurs too many times in its target genome is also removed.
     *
     * @param repId		ID of the representative genome to scan
     *
     * @return the number of hammers removed
     */
    private int cleanHammers(String repId) {
        Genome repGenome = this.repGenomes.getGenome(repId);
        log.info("Scanning genome {} for cleaning step.", repGenome);
        int retVal = 0;
        CountMap<String> counts = new CountMap<String>();
        // Scan all the contig sequences in both directions.
        for (Contig contig : repGenome.getContigs()) {
            String seq = contig.getSequence();
            String rSeq = Contig.reverse(seq);
            Set<String> hammers = this.getHammers(List.of(seq, rSeq));
            for (String kmer : hammers) {
                HammerScore score = this.hammerMap.get(kmer);
                if (score != null) {
                    // Here we have a hammer.  If it is invalid, remove it.
                    if (! repId.contentEquals(Feature.genomeOf(score.getFid()))) {
                        HammerScore testScore = this.hammerMap.remove(kmer);
                        // If the remove returned NULL, someone else got to this hammer first.
                        if (testScore != null)
                            retVal++;
                    } else {
                        // Here the hammer is in the same genome.  Count it.  Counts greater than
                        // the maximum are removed.
                        counts.count(kmer);
                    }
                }
            }
        }
        int dups = 0;
        for (var counter : counts.counts()) {
            if (counter.getCount() > this.maxCount) {
                HammerScore testScore = this.hammerMap.remove(counter.getKey());
                // If the remove returned NULL, someone else got to this hammer first.
                if (testScore != null)
                    dups++;
            }
        }
        if (dups > 0) {
            log.info("{} duplicate hammers removed from {}.", dups, repGenome);
            retVal += dups;
        }
        log.info("{} hammers cleaned using {}.", retVal, repGenome);
        return retVal;
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
     *
     * @throws IOException
     */
    private void writeHammers(PrintWriter writer) throws IOException {
        log.info("Writing output (this will take many minutes).");
        long hammersIn = 0;
        long lowPrec = 0;
        long lowWorth = 0;
        long hammersOut = 0;
        long badHammer = 0;
        // Create the scoring helper object.
        Object scoreHelper = this.strengthType.helperScan(this.hammerMap);
        // Set up the reject-file stream here.
        PrintWriter rWriter = null;
        try {
            if (this.rejectFile != null) {
                rWriter = new PrintWriter(this.rejectFile);
                log.info("Writing reject list to {}.", this.rejectFile);
                this.writeHammerHeader(rWriter);
            }
            // Write the output header.
            this.writeHammerHeader(writer);
            // Loop through the hammers in the hammer map.
            for (var hammerEntry : this.hammerMap) {
                String hammer = hammerEntry.getKey();
                HammerScore score = hammerEntry.getValue();
                // Only proceed if the hammer is good.
                hammersIn++;
                if (score.isBadHammer()) {
                    badHammer++;
                    if (rWriter != null)
                        this.writeHammer(rWriter, hammer, score, scoreHelper);
                } else {
                    double worth = score.getWorthiness();
                    double precision = score.getPrecision();
                    // Filter out the hammers that aren't useful.
                    if (worth < this.minWorth)
                        lowWorth++;
                    else if (precision < this.minPrec)
                        lowPrec++;
                    else {
                        // Here we have a good hammer.
                        this.writeHammer(writer, hammer, score, scoreHelper);
                        hammersOut++;
                    }
                }
            }
            log.info("{} hammers checked, {} marked bad, {} rejected due to low precision, {} due to low worth, {} output.",
                    hammersIn, badHammer, lowPrec, lowWorth, hammersOut);
        } finally {
            // Insure the reject file is closed.
            if (rWriter != null)
                rWriter.close();
        }
    }

    /**
     * Write a hammer to an output file.
     *
     * @param writer		output print writer
     * @param hammer		hammer to write
     * @param score			associated scoring object
     * @param scoreHelper	helper object for scoring
     */
    private void writeHammer(PrintWriter writer, String hammer, HammerScore score, Object scoreHelper) {
        writer.println(hammer + "\t" + score.getFid() + "\t" + score.getStrength(scoreHelper) + "\t"
                + score.getPrecision() + "\t" + score.getWorthiness() + "\t"
                + score.getGoodHits() + "\t" + score.getRoleId());
    }

    /**
     * Write the output header for a hammer list.
     *
     * @param writer	output writer for the list
     */
    private void writeHammerHeader(PrintWriter writer) {
        writer.println("hammer\tfid\tstrength\tprecision\tworth\thits\trole");
    }


}
