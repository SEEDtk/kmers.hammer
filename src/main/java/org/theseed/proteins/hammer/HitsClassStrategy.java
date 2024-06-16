/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Collection;

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
    public ScoreMap computeScores(Collection<HammerDb.Hit> hits, int len, double covg) {
        double weight = this.getWeight(len, covg);
        ScoreMap retVal = new ScoreMap();
        for (var hit : hits) {
            retVal.count(hit.getGenomeId(), weight * this.getWeight(hit), hit.getRole());
        }
        return retVal;
    }

}
