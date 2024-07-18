/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Set;

import org.theseed.genome.Feature;

/**
 * This hammer scoring algorithm weights a hammer by the number of hammers that yield the given role and genome ID as it compares
 * to the mean for all role / genome pairs.  Thus, if a hit result has lots of hammers, each hit has less weight than a hit result
 * with few hammers.
 *
 * @author Bruce Parrello
 *
 */
public class RoleHammerScore extends HammerScore {

    public RoleHammerScore(String fid, String role, Set<String> neighborSet, boolean badFlag) {
        super(fid, role, neighborSet, badFlag);
    }

    @Override
    public double getStrength(Object helper) {
        RoleCounters realHelper = (RoleCounters) helper;
        String genomeId = Feature.genomeOf(this.getFid());
        double retVal = realHelper.getWeight(this.getRoleId(), genomeId);
        return retVal;
    }

}
