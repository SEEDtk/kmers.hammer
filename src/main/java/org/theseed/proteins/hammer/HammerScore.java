/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Set;

import org.theseed.genome.Feature;

/**
 * This object is used to help score hammers being created using a protein finder.  It tells us
 * its source feature ID, and the good and bad hit counts.  Each hit count represents a single
 * neighbor (good) or distant (bad) genome.  Programs that use this object are required to
 * process sequences in genome order, so that all potential hammers in one genome are processed
 * before those in the next genome.  The object stores the ID of the last genome counted.  If
 * the current hit is from the same genome, the counters are not updated.
 *
 * @author Bruce Parrello
 *
 */
public class HammerScore {

    // FIELDS
    /** source feature ID */
    private String fid;
    /** ID of last genome processed */
    private String lastGenome;
    /** neighbor set for this hammer's source genome (including the source itself) */
    private Set<String> neighborhood;
    /** number of neighborhood hits */
    private int goodHits;
    /** number of foreign hits */
    private int badHits;
    /** total number of genomes in the master database */
    private static int totalGenomes = 0;
    /** mean neighborhood size */
    private static double NORMAL_NEIGHBORHOOD = 42.0;

    /**
     * Construct a new hammer scoring object.
     *
     * @param fid			source feature for the hammer
     * @param neighborSet	ID set for the source genome's neighborhood
     */
    public HammerScore(String fid, Set<String> neighborSet) {
        this.fid = fid;
        this.neighborhood = neighborSet;
        this.lastGenome = "";
        // We start with 1 good hit from the originating representative genome.  There are no bad hits.
        // A hammer with a bad hit in the repgen set is deleted from the potential-hammer set immediately.
        this.goodHits = 1;
        this.badHits = 0;
    }

    /**
     * Record a collision with this hammer.
     *
     * @param fid	ID of the colliding feature's genome
     */
    public void recordHit(String genomeId) {
        if (! genomeId.contentEquals(this.lastGenome)) {
            // Here we need to count the hit.
            if (this.neighborhood.contains(genomeId))
                this.goodHits++;
            else
                this.badHits++;
            // Insure this genome isn't counted again.
            this.lastGenome = genomeId;
        }
    }

    /**
     * @return TRUE if the specified feature is in a different genome from this hammer, else FALSE
     *
     * @param fid2	ID of hammer's alternative source feature
     */
    public boolean isDisqualifyingHit(String fid2) {
        boolean retVal = ! (Feature.genomeOf(this.fid).contentEquals(Feature.genomeOf(fid2)));
        return retVal;
    }

    /**
     * Specify the total number of genomes being scanned.  This is
     * determined once at the beginning of the run.
     *
     * @param total		total number of genomes being scanned
     */
    public static void setTotalGenomes(int total) {
        totalGenomes = total;
    }

    /**
     * @return the precision of this hammer
     */
    public double getPrecision() {
        double retVal;
        if (this.badHits == 0)
            retVal = 1.0;
        else {
            int distantGenomes = totalGenomes - this.neighborhood.size();
            retVal = (distantGenomes - this.badHits) / (double) distantGenomes;
        }
        return retVal;
    }

    /**
     * @return the worthiness of this hammer
     */
    public double getWorthiness() {
        return this.goodHits / (double) this.neighborhood.size();
    }

    /**
     * The strength is a measure of our confidence in the hammer.  We start with a fixed value of 1/42
     * (42 being the mean neighborhood size), then add the fraction of the neighborhood outside the
     * representative genome hit by the hammer.  The idea is to give extra weight to evidence that the
     * hammer has worth, without dropping singleton neighborhoods to a weight of 0.
     *
     * @return the strength of the hammer
     */
    public double getStrength() {
        double retVal = 1.0 / NORMAL_NEIGHBORHOOD;
        if (this.goodHits > 1)
            retVal += (this.goodHits - 1);
        return retVal / this.neighborhood.size();
    }

    /**
     * @return the hammer's source feature ID
     */
    public String getFid() {
        return this.fid;
    }

    /**
     * @return the number of good hits for this hammer
     */
    public int getGoodHits() {
        return this.goodHits;
    }

}
