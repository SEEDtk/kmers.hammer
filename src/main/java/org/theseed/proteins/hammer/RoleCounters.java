/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.HashMap;
import java.util.Map;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;

/**
 * This object is the score helper for role-based scoring.  It scans the scores in the hammer map and counts the
 * number of hammers for each role/genome pair.  It can return the mean for all pairs and the count for an individual
 * pair.  The mean will be nonzero if the hammer map is non-empty.
 *
 * We set this up with a count map for each role, keyed by genome ID.
 *
 *
 * @author Bruce Parrello
 *
 */
public class RoleCounters {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(RoleCounters.class);
    /** count hash */
    private Map<String, CountMap<String>> roleCountMaps;
    /** mean count value */
    private double meanCount;

    /**
     * Scan the hammer map to produce the counts and compute the mean.
     *
     * @param hammers	hammer map to scan
     */
    public RoleCounters(HammerMap<HammerScore> hammers) {
        log.info("Scanning hammer map to compute role counts.");
        long lastMsg = System.currentTimeMillis();
        this.roleCountMaps = new HashMap<String, CountMap<String>>(40);
        long count = 0;
        for (var hammerEntry : hammers) {
            count++;
            HammerScore score = hammerEntry.getValue();
            String genomeId = Feature.genomeOf(score.getFid());
            String role = score.getRoleId();
            CountMap<String> countMap = this.roleCountMaps.computeIfAbsent(role, x -> new CountMap<String>());
            countMap.count(genomeId);
            long now = System.currentTimeMillis();
            if (now - lastMsg >= 5000) {
                log.info("{} hammers scanned for role counting.", count);
                lastMsg = now;
            }
        }
        log.info("{} total hammers scanned for role counting.", count);
        // Now compute the mean.
        count = 0;
        long total = 0;
        for (CountMap<String> map : this.roleCountMaps.values()) {
            for (var counter : map.counts()) {
                total += counter.getCount();
                count++;
            }
        }
        if (total == 0)
            this.meanCount = 1.0;
        else
            this.meanCount = ((double) total) / count;
        log.info("Mean count is {}.", this.meanCount);
    }

    /**
     * @return the weighted score for a given role / genome pair
     *
     * @param 	roleId		ID of the relevant role
     * @param	genomeId	ID of the relevant genome
     */
    public double getWeight(String roleId, String genomeId) {
        int pairCount;
        CountMap<String> countMap = this.roleCountMaps.get(roleId);
        if (countMap == null)
            pairCount = 0;
        else
            pairCount = countMap.getCount(genomeId);
        double retVal;
        if (pairCount == 0)
            retVal = 0.0;
        else
            retVal = this.meanCount / pairCount;
        return retVal;
    }

}
