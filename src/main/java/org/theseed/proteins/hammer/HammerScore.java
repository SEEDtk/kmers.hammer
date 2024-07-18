/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Set;

/**
 * This object is used to help score hammers being created using a protein finder.  It tells us
 * its source feature ID, and the good and bad hit counts.  Each hit count represents a single
 * neighbor (good) or distant (bad) genome.  Each hammer is only processed once per genome per
 * role.  If a hammer is found in multiple roles, it is marked as bad.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerScore implements HammerMap.IScore {

    // FIELDS
    /** source feature ID */
    private String fid;
    /** role for this hammer */
    private String roleId;
    /** neighbor set for this hammer's source genome (including the source itself) */
    protected Set<String> neighborhood;
    /** number of neighborhood hits */
    protected int goodHits;
    /** number of foreign hits */
    protected int badHits;
    /** TRUE if the hammer is invalid (found in multiple repgens), else FALSE */
    private boolean badHammer;
    /** mean neighborhood size */
    private static double NORMAL_NEIGHBORHOOD = 42.0;
    /** total number of genomes in the system */
    protected static int totalGenomes;

    /**
     * This enum handles the different types of scoring objects.  Each uses a different method for computing
     * the strength of a hammer (which is our measure of confidence.
     *
     * @author Bruce Parrello
     *
     */
    public enum Type {

        /** strength is related to worthiness; a hammer that occurs in more neighbors is stronger */
        WORTH_BASED {
            @Override
            public HammerScore create(String fid, String role, Set<String> neighborSet, boolean badFlag) {
                return new HammerScore.Worth(fid, role, neighborSet, badFlag);
            }

            @Override
            public Object helperScan(HammerMap<HammerScore> hammers) {
                return null;
            }
        },

        /** strength is related to neighbor size; a hammer for a more-common genome is stronger */
        POP_BASED {
            @Override
            public HammerScore create(String fid, String role, Set<String> neighborSet, boolean badFlag) {
                return new HammerScore.Pop(fid, role, neighborSet, badFlag);
            }

            @Override
            public Object helperScan(HammerMap<HammerScore> hammers) {
                return null;
            }
        },

        /** strength is determined by the ratio between the number of hammers of a given type and the mean */
        RATIO_BASED {
            @Override
            public HammerScore create(String fid, String role, Set<String> neighborSet, boolean badFlag) {
                return new RoleHammerScore(fid, role, neighborSet, badFlag);
            }

            @Override
            public Object helperScan(HammerMap<HammerScore> hammers) {
                return new RoleCounters(hammers);
            }
        };

        /**
         * This creates a scoring object for a hammer.
         *
         * @param fid			source feature for the hammer
         * @param role			ID of the role for the hammer
         * @param neighborSet	ID set for the source genome's neighborhood
         * @param badFlag		initial value for the bad-hammer flag
         *
         * @return a hammer-scoring object of this type
         */
        public abstract HammerScore create(String fid, String role, Set<String> neighborSet, boolean badFlag);

        /**
         * This scans the hammer map to produce a helper object for scoring.
         *
         * @param hMap		hammer map to scan
         */
        public abstract Object helperScan(HammerMap<HammerScore> hammers);

    }

    /**
     * This hammer type has more confidence if the genome type is more common.
     */
    public static class Pop extends HammerScore {

        public Pop(String fid, String role, Set<String> neighborSet, boolean badFlag) {
            super(fid, role, neighborSet, badFlag);
        }

        @Override
        public double getStrength(Object helper) {
            double retVal = this.neighborhood.size() / NORMAL_NEIGHBORHOOD;
            if (retVal > 1.0) retVal = 1.0;
            return retVal;
        }

    }

    /**
     * This hammer type has more confidence if the hammer itself is more common.
     */
    public static class Worth extends HammerScore {

        public Worth(String fid, String role, Set<String> neighborSet, boolean badFlag) {
            super(fid, role, neighborSet, badFlag);
        }

        @Override
        public double getStrength(Object helper) {
            double retVal = 1.0 / NORMAL_NEIGHBORHOOD;
            if (this.goodHits > 1)
                retVal += (this.goodHits - 1);
            return retVal / this.neighborhood.size();
        }

    }


    /**
     * Construct a new hammer scoring object.
     *
     * @param fid			source feature for the hammer
     * @param role			ID of the role for the hammer
     * @param neighborSet	ID set for the source genome's neighborhood
     * @param badFlag		initial value for the bad-hammer flag
     */
    public HammerScore(String fid, String role, Set<String> neighborSet, boolean badFlag) {
        this.fid = fid;
        this.neighborhood = neighborSet;
        this.badHammer = badFlag;
        this.roleId = role;
        // There is one good hit for the representative genome where we found the hammer, and no bad hits.
        this.goodHits = 1;
        this.badHits = 0;
    }

    /**
     * Specify the total number of genomes in the system.
     *
     * @param total		value to specify
     */
    public static void setTotalGenomes(int total) {
        totalGenomes = total;
    }

    /**
     * Record a hit with this hammer.
     *
     * @param fid		ID of the colliding feature's genome
     * @param role		ID of the role in which the hammer was found
     */
    public synchronized void recordHit(String genomeId, String role) {
        if (! this.roleId.contentEquals(role))
            this.badHammer = true;
        else if (this.neighborhood.contains(genomeId))
            this.goodHits++;
        else
            this.badHits++;
    }

    /**
     * @return the precision of this hammer
     */
    public double getPrecision() {
        double retVal;
        if (this.badHits == 0)
            retVal = 1.0;
        else {
            int distantGenomes = (totalGenomes - this.neighborhood.size());
            retVal = (distantGenomes - this.badHits) / (double) (distantGenomes);
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
     * The strength is a measure of our confidence in the hammer.
     *
     * @param helper	helper object produced by the hammer scan
     *
     * @return the strength of the hammer
     */
    public abstract double getStrength(Object helper);

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

    /**
     * @return TRUE if this hammer fails one of the validity tests
     */
    @Override
    public boolean isBadHammer() {
        return this.badHammer;
    }

    /**
     * Soecify that this hammer is bad.
     */
    @Override
    public void setBadHammer() {
        this.badHammer = true;
    }

    /**
     * @return the ID of the hammer's role
     */
    public String getRoleId() {
        return this.roleId;
    }


}
