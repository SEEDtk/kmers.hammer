/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.locations.LocationFinder;
import org.theseed.reports.NaturalSort;
import org.theseed.sequence.RnaKmers;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command scans the results of a synthetic-sample hammer test and analyzes individual bad hits.  The command line
 * specifies the genome ID whose hits are to be analyzed.  For each hit, we determine the identity of the hammer, plus
 * the ID and role of the source feature for the hammer and the target feature hit (if any).
 *
 * The standard input should contain the result file from the contig test.  The report will be written to the standard
 * output.  The positional parameters are the ID of the genome of interest and the name of the genome source (file or
 * directory) containing all of the genomes processed.  The report will only include hits against the specified genome.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	name of input file (if not STDIN)
 * -o	name of output file (if not STDOUT)
 *
 * --source		type of genome source (default DIR)
 *
 *
 * @author Bruce Parrello
 *
 */
public class ContigTestLocationProcessor extends BasePipeProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ContigTestLocationProcessor.class);
    /** genome source */
    private GenomeSource genomes;
    /** target genome */
    private Genome targetGenome;
    /** location finder for target genome */
    private LocationFinder targetFinder;
    /** input column index for the hammer hit location */
    private int hitColIdx;
    /** input column index for the hammer feature ID */
    private int fidColIdx;
    /** input column index for the expected representative genome ID */
    private int repColIdx;
    /** current hammer genome */
    private Genome hammerGenome;
    /** expected representative genome */
    private Genome repGenome;
    /** number of hits processed */
    private int hitCount;
    /** position in a hammer hit location input string at which to find the real location string */
    private int residualIdx;
    /** set of hits being kept */
    private Set<HammerHit> hitSet;
    /** comparator for natural ordering */
    private static final Comparator<String> SORTER = new NaturalSort();

    // COMMAND-LINE OPTIONS

    /** type of genome source (MASTER, CORE, DIR, etc.) */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** ID of the genome of interest */
    @Argument(index = 0, metaVar = "genomeID", usage = "ID of the genome whose hammer hits are to be analyzed", required = true)
    private String targetID;

    /** name of the genome source file or directory */
    @Argument(index = 1, metaVar = "sourceDir", usage = "name of the file or directory containing the genomes needed", required = true)
    private File sourceDir;

    /**
     * This class tracks data on a single hit.  It is sorted by source hammer genome and feature ID.
     */
    private class HammerHit implements Comparable<HammerHit> {

        /** location of target genome hit */
        private Location hitLoc;
        /** sequence of the hammer */
        private String hammerDna;
        /** feature hit by the hammer */
        private Set<Feature> hitFeatures;
        /** genome ID of hammer source */
        private String sourceGenomeId;
        /** feature ID of hammer source */
        private String sourceFeatureId;
        /** TRUE for a good hit, FALSE for a bad hit */
        private boolean goodHit;
        /** index number of this hit, to guarantee each hit is unique */
        private int idx;

        /**
         * Construct a hammer hit from input data.
         *
         * @param hitLocString	hit location string (excluding genome ID prefix)
         * @param hammerFid		feature ID of the hammer
         * @param repId			expected representative genome ID
         */
        protected HammerHit(String hitLocString, String hammerFid, String repId) {
            // Parse the hit location.
            this.hitLoc = Location.parseSeedLocation(hitLocString);
            this.hammerDna = ContigTestLocationProcessor.this.targetGenome.getDna(this.hitLoc);
            // Compute the features hit at the location.
            this.hitFeatures = ContigTestLocationProcessor.this.targetFinder.getFeatures(this.hitLoc);
            // Extract the hammer source genome ID.
            this.sourceGenomeId = Feature.genomeOf(hammerFid);
            this.sourceFeatureId = hammerFid;
            // Determine if this is a good hit.
            this.goodHit = (this.sourceGenomeId.contentEquals(repId));
            // Save the index number.
            this.idx = ContigTestLocationProcessor.this.hitCount;
            ContigTestLocationProcessor.this.hitCount++;
        }


        @Override
        public int compareTo(HammerHit o) {
            // Put bad hits first.
            int retVal = Boolean.compare(this.goodHit, o.goodHit);
            if (retVal == 0) {
                // The natural sort will sort by taxon, then genome version number, then peg number.
                retVal = SORTER.compare(this.sourceFeatureId, o.sourceFeatureId);
                // Fall back to sort by hit location.
                if (retVal == 0)
                    retVal = this.idx - o.idx;
            }
            return retVal;
        }

        /**
         * @return the DNA of the hammer that hit the target
         */
        protected String getHammerDna() {
            return this.hammerDna;
        }

        /**
         * @return the set of features at the hit location
         */
        protected Set<Feature> getHitFeatures() {
            return this.hitFeatures;
        }

        /**
         * @return the ID of the genome containing the original hammer
         */
        protected String getSourceGenomeId() {
            return this.sourceGenomeId;
        }

        /**
         * @return the ID of the feature containing the original hammer
         */
        protected String getSourceFeatureId() {
            return this.sourceFeatureId;
        }

        /**
         * @return TRUE if the hit found the correct source genome
         */
        protected boolean isGoodHit() {
            return this.goodHit;
        }

        /**
         * @return the hammer length
         */
        protected int getHammerLen() {
            return this.hitLoc.getLength();
        }

    }

    @Override
    protected void setPipeDefaults() {
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Set up the genome source.
        if (! this.sourceDir.exists())
            throw new FileNotFoundException("Genome source " + this.sourceDir + " does not exist.");
        log.info("Connecting to genome source {} of type {}.", this.sourceDir, this.sourceType);
        this.genomes = this.sourceType.create(this.sourceDir);
        if (! this.genomes.getIDs().contains(this.targetID))
            throw new ParseFailureException("Target genome " + this.targetID + " not found in genome source " + this.sourceDir);
        // Load the target genome.
        log.info("Loading target genome from {}.", this.sourceDir);
        this.targetGenome = this.genomes.getGenome(this.targetID);
        this.targetFinder = new LocationFinder(this.targetGenome);
        // Denote we have not loaded the representative genome.
        this.repGenome = null;
        // Compute the position in each hit location string from the input where the real location string begins.
        this.residualIdx = this.targetID.length() + 1;
        // Initialize the hit tracking.
        this.hitCount = 0;
        this.hitSet = new TreeSet<HammerHit>();
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Verify that we have the input columns we need.
        this.hitColIdx = inputStream.findField("location");
        this.fidColIdx = inputStream.findField("hammer_fid");
        this.repColIdx = inputStream.findField("rep_id");
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // We loop through the input file, collecting hits to our genome.  These are sorted by source
        // hammer, so that when we process them the output is more organized.
        log.info("Reading input stream.");
        int inCount = 0;
        for (var line : inputStream) {
            // Get the location of the hammer hit.
            String hammerHitLocString = line.get(this.hitColIdx);
            // Split out the genome ID.
            String hammerGenomeId = StringUtils.substringBefore(hammerHitLocString, ":");
            if (hammerGenomeId.contentEquals(this.targetID)) {
                // Here we are interested in this hit.
                String hitLocString = hammerHitLocString.substring(this.residualIdx);
                String hammerFid = line.get(this.fidColIdx);
                String repId = line.get(this.repColIdx);
                HammerHit hit = this.new HammerHit(hitLocString, hammerFid, repId);
                this.hitSet.add(hit);
                // Insure the representative genome is in memory.
                if (this.repGenome == null) {
                    this.repGenome = this.genomes.getGenome(repId);
                    log.info("Representative genome {} loaded.", this.repGenome);
                }
            }
            // Count the line in.
            inCount++;
            if (log.isInfoEnabled() && inCount % 10000 == 0)
                log.info("{} lines read, {} hits kept.", inCount, this.hitCount);
        }
        log.info("{} total hammer hits kept.", this.hitCount);
        // Now we process the individual hits to produce the output.
        log.info("Writing output.");
        writer.println("target_fid\ttarget_function\tsim_to_hammer\tsim_to_rep\thammer\thammer_genome\thammer_fid\thammer_function\tgood_hit");
        // These will track the current hammer feature and genome.
        String hammerFid = "";
        String hammerGenomeId = "";
        Feature hammerFeature = null;
        String hammerFunction = "";
        RnaKmers repKmers = null;
        RnaKmers hammerKmers = null;
        // We need some counters.
        int goodHits = 0;
        int badHits = 0;
        int multiHits = 0;
        int openHits = 0;
        for (HammerHit hit : this.hitSet) {
            // Remember the hammer length.
            int hammerLen = hit.getHammerLen();
            // Compute the hammer's source feature.  This may require loading a new hammer genome.
            if (! hammerGenomeId.contentEquals(hit.getSourceGenomeId())) {
                hammerGenomeId = hit.getSourceGenomeId();
                this.hammerGenome = this.genomes.getGenome(hammerGenomeId);
                log.info("Hammer genome {} loaded.", this.hammerGenome);
            }
            if (! hammerFid.contentEquals(hit.getSourceFeatureId())) {
                // Here we have a new hammer feature.  This is a feature in the genome the hammer came from.
                hammerFid = hit.getSourceFeatureId();
                hammerFeature = this.hammerGenome.getFeature(hammerFid);
                hammerFunction = hammerFeature.getPegFunction();
                // Get a kmer structure for the corresponding features in the representative genome.  We scan the whole
                // genome and pick off each sequence that has the same function.
                repKmers = new RnaKmers(hammerLen);
                for (Feature repFeat : this.repGenome.getPegs()) {
                    if (repFeat.getPegFunction().contentEquals(hammerFunction))
                        repKmers.addSequence(repFeat.getDna());
                }
                // Get a kmer structure for the hammer feature.
                hammerKmers = new RnaKmers(hammerFeature.getDna(), hammerLen);
            }
            // Get the hammer itself.
            String hammerDNA = hit.getHammerDna();
            // Determine the type of hit.
            boolean goodHit = hit.isGoodHit();
            if (goodHit)
                goodHits++;
            else
                badHits++;
            // Form the string for the trailing columns (which are constant for all the output lines of this hit).
            String trailer = "\t" + hammerDNA + "\t" + hammerGenomeId + "\t" + hit.getSourceFeatureId() + "\t"
                    + hammerFunction + "\t" + (hit.isGoodHit() ? "Y" : "");
            // Get the features hit.
            var targetFeats = hit.getHitFeatures();
            if (targetFeats.size() == 0) {
                // Here the target region is outside any called features.
                writer.println(this.targetID + "\t<uncalled region>\t\t" + trailer);
                openHits++;
            } else {
                if (targetFeats.size() > 1)
                    multiHits++;
                // Here we have one or more features covering the region hit.
                for (Feature feat : targetFeats) {
                    // Compute the sims for this feature.
                    RnaKmers featKmers = new RnaKmers(feat.getDna(), hammerLen);
                    int simToHammer = featKmers.similarity(hammerKmers);
                    int simToRep = featKmers.similarity(repKmers);
                    writer.format("%s\t%s\t%d\t%d%s%n", feat.getId(), feat.getFunction(), simToHammer, simToRep, trailer);
                }
            }
        }
        log.info("{} good hits, {} bad hits, {} open-region hits, {} multi-feature hits.", goodHits, badHits, openHits, multiHits);
    }



}
