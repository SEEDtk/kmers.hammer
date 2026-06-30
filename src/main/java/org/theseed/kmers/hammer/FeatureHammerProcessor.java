package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.proteins.hammer.HammerMap;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;

/**
 * This command analyzes potential feature hammers for reasons they were rejected.  It can be used to determine why a particular
 * feature has more or fewer hammers than expected. The analysis will be written to the standard output.
 * 
 * The positional parameters are the name of the finder directory used to generate the hammers, the ID of the feature to analyze,
 * the name of the feature's protein, and the name of the hammer file containing the hammers.
 * 
 * The command-line options are as follows.
 * 
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input file containing the repgen list data (if not STDIN)
 * -o	output file for the hammer list (if not STDOUT)
 * -K	kmer size for a hammer (default 20)
 * -t	type of genome source for the optional cleaning step (default DIR)
 *
 * --maxRepeat	maximum percent of a hammer that can belong to a single base pair type (default 0.70)
 * --clean      the directory containing the genomes to use for cleaning the hammers (default none)
 * 
 */
public class FeatureHammerProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    private static final Logger log = LoggerFactory.getLogger(FeatureHammerProcessor.class);
    /** protein finder */
    private ProteinFinder finder;
    /** cleaning source */
    private GenomeSource repGenomes;
    /** name of FASTA file for the protein */
    private File proteinDnaFasta;
    /** set of hammers accepted for the protein */
    private Set<String> acceptedHammers;
    /** record of rejected hammers by rejection reason */
    private CountMap<String> reasonReport;

   // COMMAND-LINE OPTIONS

    /** hammer size, in base pairs */
    @Option(name = "--kSize", aliases = { "-K", "--kmer" }, metaVar = "21", usage = "hammer DNA kmer size")
    private int kmerSize;

    /** maximum fraction of a hammer allowed for a particular nucleotide */
    @Option(name = "--maxRepeat", metaVar = "0.60", usage = "maximum fraction of a hammer allowed for any one nucleotide type")
    private double maxRepeat;

    /** type of genome source for the optional cleaning directory */
    @Option(name = "--type", aliases = { "-t" }, usage = "type of genome source for optional cleaning genomes")
    private GenomeSource.Type sourceType;

    /** repgen genome directory */
    @Option(name = "--clean", metaVar = "repDir", usage = "directory containing genomes to use for cleaning the hammers")
    private File repDir;

    /** protein-finder directory */
    @Argument(index = 0, metaVar = "finderDir", usage = "protein finder directory", required = true)
    private File finderDir;

    /** feature ID to analyze */
    @Argument(index = 1, metaVar = "featureID", usage = "feature ID to analyze", required = true)
    private String targetFid;

    /** protein name for the feature */
    @Argument(index = 2, metaVar = "proteinName", usage = "protein name for the feature", required = true)
    private String proteinName;

    /** hammer file */
    @Argument(index = 3, metaVar = "hammerFile", usage = "file containing the relevant hammers", required = true)
    private File hammerFile;

    @Override
    protected void setReporterDefaults() {
        this.kmerSize = 20;
        this.maxRepeat = 0.70;
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected void validateReporterParms() throws ParseFailureException, IOException {
        // Check the tuning parameters.
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer size must be at least 2.");
        if (this.maxRepeat <= 0.25 || this.maxRepeat > 1.0)
            throw new ParseFailureException("Maximum repeat fraction must be between 0.25 and 1.0.");
        // Set up the cleaning directory.
        if (this.repDir == null) {
            log.info("No cleaning step will be performed.");
            this.repGenomes = null;
        } else {
            log.info("Hammers will be cleaned using representative genomes in {}.", this.repDir);
            this.repGenomes = this.sourceType.create(this.repDir);
        }
        // Save the kmer size for the DNA and Genome kmer classes.
        DnaKmers.setKmerSize(this.kmerSize);
        GenomeKmers.setKmerSize(this.kmerSize);
        // Set up the protein finder.
        if (! this.finderDir.isDirectory())
            throw new FileNotFoundException("Finder directory " + this.finderDir + " is not found or invalid.");
        log.info("Loading protein-finder from {}.", this.finderDir);
        this.finder = new ProteinFinder(this.finderDir);
        // Get the map of protein names to DNA FASTA files.
        Map<String, File> roleMap = this.finder.getFastas();
        // Save the FASTA file for the named protein.
        this.proteinDnaFasta = roleMap.get(this.proteinName);
        if (this.proteinDnaFasta == null)
            throw new ParseFailureException("Protein " + this.proteinName + " not found in finder directory " + this.finderDir);
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Initialize the rejection-reason count.
        this.reasonReport = new CountMap<>();
        // The first task is to run through the hammer file and extract the hammers for the specified feature. We save these in memory
        // so we know which ones were kept.
        this.acceptedHammers = this.findActualHammers();
        if (this.acceptedHammers.isEmpty())
            throw new IOException("No hammers were found for feature " + this.targetFid + " in hammer file " + this.hammerFile);
        else
            log.info("{} hammers were found for feature {}.", this.acceptedHammers.size(), this.targetFid);
        // Start the output report.
        writer.println("hammer\tcategory\tdetails");
        // We now need to perform a bunch of processing using the finder file. We will load the finder DNA into memory
        // while doing this, so a subroutine is used to insure the DNA memory is released when we are done.
        Set<String> potentialHammers = this.findPotentialHammers(writer);
        // Now we loop through the remaining potential hammers and check them for lint (the repeat fraction).
        int rejectCount = 0;
        int maxBaseCount = (int) (this.maxRepeat * this.kmerSize);
        log.info("Checking {} potential hammers for common base occurring more than {} times.", potentialHammers.size(), maxBaseCount);
        Iterator<String> iter = potentialHammers.iterator();
        while (iter.hasNext()) {
            String hammer = iter.next();
            int bestBaseCount = HammerMap.commonBaseCount(hammer);
            if (bestBaseCount > maxBaseCount) {
                this.reportReason(writer, hammer, "LINT", "most common base occurs " + bestBaseCount + " times");
                iter.remove();
                rejectCount++;
            }
        }
        log.info("Rejected {} hammers for lint.", rejectCount);
        // Only proceed if we have potential hammers left and we have a cleaning genome source.
        if (this.repGenomes != null && ! potentialHammers.isEmpty()) {
            // Now we check the remaining potential hammers against the representative genomes. We do one pass to get
            // everything but the source genome, then do the source genome last, since presence in the source genome
            // has a different reason code.
            String targetGenomeId = Feature.genomeOf(this.targetFid);
            Set<String> genomeIds = this.repGenomes.getIDs();
            rejectCount = 0;
            for (String genomeId : genomeIds) {
                if (! genomeId.contentEquals(targetGenomeId)) {
                    log.info("Checking {} potential hammers against genome {}.", potentialHammers.size(), genomeId);
                    GenomeKmers genomeKmers = new GenomeKmers(this.repGenomes.getGenome(genomeId), this.kmerSize);
                    iter = potentialHammers.iterator();
                    while (iter.hasNext()) {
                        String hammer = iter.next();
                        if (genomeKmers.contains(hammer)) {
                            // This hammer is found in the genome.  We will reject it.
                            iter.remove();
                            this.reportReason(writer, hammer, "CLEAN", "found in cleaning genome " + genomeId);
                            rejectCount++;
                        }
                    }
                }
            }
            log.info("Rejected {} potential hammers for being found in cleaning genomes.", rejectCount);
            // Now we check the source genome.
            Genome sourceGenome = this.repGenomes.getGenome(targetGenomeId);
            if (sourceGenome == null)
                log.error("Source genome is not found in repgen set!");
            else {
                log.info("Checking {} potential hammers against source genome {}.", potentialHammers.size(), targetGenomeId);
                // Now we look at the source genome and pull out the DNA sequences with the target feature removed.
                Feature feat = sourceGenome.getFeature(this.targetFid);
                if (feat == null)
                    throw new IOException("Source genome " + targetGenomeId + " does not contain feature " + this.targetFid + ".");
                Location featLoc = feat.getLocation();
                String featContigId = featLoc.getContigId();
                // We create a list of the contig sequences. When we come to the target feature's contig, we output the sequence on either side.
                List<DnaKmers> contigKmers = new ArrayList<>();
                for (Contig contig : sourceGenome.getContigs()) {
                    String seq = contig.getSequence();
                    if (contig.getId().contentEquals(featContigId)) {
                        // This is the target feature's contig.  We need to remove the feature from it.
                        int start = featLoc.getLeft();
                        int end = featLoc.getRight();
                        String leftSeq = seq.substring(0, start);
                        String rightSeq = seq.substring(end);
                        if (! leftSeq.isBlank())
                            contigKmers.add(new DnaKmers(leftSeq, this.kmerSize));
                        if (! rightSeq.isBlank())
                            contigKmers.add(new DnaKmers(rightSeq, this.kmerSize));
                    } else {
                        // This is a normal contig.  We just add it.
                        contigKmers.add(new DnaKmers(seq, this.kmerSize));
                    }
                }
                rejectCount = 0;
                iter = potentialHammers.iterator();
                while (iter.hasNext()) {
                    String hammer = iter.next();
                    // Check this hammer against each contig.
                    boolean found = false;
                    Iterator<DnaKmers> contigIter = contigKmers.iterator();
                    while (! found && contigIter.hasNext()) {
                        DnaKmers contigKmer = contigIter.next();
                        if (contigKmer.contains(hammer))
                            found = true;
                    }
                    if (found) {
                        // This hammer is found in the genome.  We will reject it.
                        iter.remove();
                        this.reportReason(writer, hammer, "SOURCE", "found in source genome");
                        rejectCount++;
                    }
                }
                log.info("Rejected {} potential hammers for being found in source genome.", rejectCount);
            }
        }
        // Finally, we reject the remaining potential hammers for reason "OTHER". Usually, this means they failed the anchor test.
        for (String hammer : potentialHammers) {
            this.reportReason(writer, hammer, "OTHER", "miscellaneous reason");
        }
        // Now we summarize the rejection reasons.
        log.info("Rejection summary:");
        for (var reasonCount : this.reasonReport.sortedCounts()) {
            log.info("{}: {} hammers", reasonCount.getKey(), reasonCount.getCount());
        }
    }

    /**
     * This method reads the hammer file and extracts the hammers that were accepted for the current feature.
     * 
     * @return a set of the accepted hammers
     * 
     * @throws IOException
     */
    private Set<String> findActualHammers() throws IOException {
        Set<String> retVal = new HashSet<>();
        // Open the hammer file. The feature ID is in column "fid", and the hammer text in column "hammer".
        try (TabbedLineReader hammerStream = new TabbedLineReader(this.hammerFile)) {
            int fidCol = hammerStream.findField("fid");
            int hammerCol = hammerStream.findField("hammer");
            // We set up some counters to show progress.
            int inCount = 0;
            int keptCount = 0;
            log.info("Reading hammer file {}.", this.hammerFile);
            long lastMessage = System.currentTimeMillis();
            // Read the hammer file line by line, saving the hammers from the specified feature ID.
            for (TabbedLineReader.Line line : hammerStream) {
                inCount++;
                if (line.get(fidCol).contentEquals(this.targetFid)) {
                    retVal.add(line.get(hammerCol));
                    keptCount++;
                }
                long now = System.currentTimeMillis();
                if (now - lastMessage > 5000) {
                    log.info("Read {} lines, kept {} hammers.", inCount, keptCount);
                    lastMessage = now;
                }
            }
        }
        return retVal;
    }

    /**
     * This method loads the DNA for the various proteins into memory and then finds the potential hammers for the specified feature.
     * Potentials that conflict with other sequences in the finder file are eliminated by this method and detailed in the output
     * report.
     * 
     * @param writer    output print writer for reporting rejected hammers
     * 
     * @return a set of the potential hammers
     * 
     * @throws IOException
     */
    private Set<String> findPotentialHammers(PrintWriter writer) throws IOException {
        // We will store the DNA sequences in here.
        Map<String, String> dnaMap = new HashMap<>();
        // We will save the target feature DNA sequence in here.
        String targetDna = null;
        // Open the DNA FASTA file and read the sequences into memory.
        try (FastaInputStream fastaStream = new FastaInputStream(this.proteinDnaFasta)) {
            log.info("Reading DNA FASTA file {}.", this.proteinDnaFasta);
            int count = 0;
            for (Sequence record : fastaStream) {
                String fid = record.getLabel();
                String dna = record.getSequence();
                count++;
                if (fid.contentEquals(this.targetFid))
                    targetDna = dna;
                else
                    dnaMap.put(fid, dna);
            }
            log.info("Read {} DNA sequences from FASTA file.", count);
            if (targetDna == null)
                throw new IOException("Target feature " + this.targetFid + " not found in FASTA file " + this.proteinDnaFasta);
        }
        // Now we get the initial set of potential hammers for the target feature.
        DnaKmers targetKmers = new DnaKmers(targetDna, this.kmerSize);
        Set<String> retVal = targetKmers.getKmerSet();
        // Remove the accepted hammers from the potential set.  These are already known to be good.
        retVal.removeAll(this.acceptedHammers);
        // The next step is to run through the DNA sequences for the other features and eliminate any un-accepted hammers that are
        // found in them. We will also report the reason for the rejection.
        int rejectCount = 0;
        for (Map.Entry<String, String> entry : dnaMap.entrySet()) {
            String fid = entry.getKey();
            String dna = entry.getValue();
            DnaKmers otherKmers = new DnaKmers(dna, this.kmerSize);
            Iterator<String> iter = retVal.iterator();
            while (iter.hasNext()) {
                String hammer = iter.next();
                if (otherKmers.contains(hammer)) {
                    // This hammer is found in another feature.  We will reject it.
                    iter.remove();
                    this.reportReason(writer, hammer, "FINDER", "found in other " + this.proteinName + " feature " + fid);
                    rejectCount++;
                }
            }
        }
        log.info("{} potential hammers were rejected because they were found in other features. {} remaining.",
            rejectCount, retVal.size());
        return retVal;
    }

    /**
     * Record the reason for a hammer rejection. We write it to the output report and also count it for a 
     * summary at the end.
     * 
     * @param writer    output report writer
     * @param hammer    hammer being rejected
     * @param code      reason code for the rejection
     * @param details    additional details for the rejection
     */
    private void reportReason(PrintWriter writer, String hammer, String code, String details) {
        writer.println(hammer + "\t" + code + "\t" + details);
        this.reasonReport.count(code);
    }

}
