/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Collection;

import org.theseed.counters.WeightMap;

/**
 * This is the simplest classification strategy for a hammer sequence.  We simply count the hits for each genome.
 * The hits are weighted by sequence length and coverage, scaled by the typical read length.
 *
 * @author Bruce Parrello
 *
 */
public class HitsClassStrategy extends ClassStrategy {

    public HitsClassStrategy(IParms processor) {
        super(processor);
    }

    @Override
    public WeightMap computeScores(Collection<HammerDb.Hit> hits, int len, double covg) {
        double weight = this.getWeight(len, covg);
        WeightMap retVal = new WeightMap();
        for (var hit : hits)
            retVal.count(hit.getGenomeId(), weight);
        return retVal;
    }

}
