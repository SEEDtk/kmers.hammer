/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Collection;

import org.theseed.proteins.hammer.HammerDb.Hit;

/**
 * This method assigns only a single classification group to a sequence.  The group with the majority of hits
 * is given a score weighted by the length and coverage (scaled by the typical length of a read.
 *
 * @author Bruce Parrello
 *
 */
public class RegionsClassStrategy extends ClassStrategy {

    // FIELDS
    /** minimum additional score required to qualify as the best group */
    private double minScore;
    /** hammer counting method */
    private HammerDb.Method method;

    public RegionsClassStrategy(IParms processor) {
        super(processor);
        // Save the threshold score and the count method.
        this.minScore = processor.getMinDiff();
        this.method = processor.getMethod();

    }

    @Override
    public ScoreMap computeScores(Collection<Hit> hits, int len, double covg) {
        // This will be the return map.
        ScoreMap retVal = new ScoreMap(2);
        // Count the genomes.
        var counts = new ScoreMap();
        for (var hit : hits)
            counts.count(hit.getGenomeId(), method.getWeight(hit), hit.getRole());
        // Only continue if something was hit.
        if (counts.size() > 0) {
            // Get the sorted counts.
            var sortedCounts = counts.sortedCounts();
            // Extract the top two counts.  If there is only one, the second count is 0.
            var bestCounter = sortedCounts.get(0);
            var secondCount = (sortedCounts.size() <= 1 ? 0.0 : sortedCounts.get(1).getCount());
            // If the threshold is met, count the best group.
            if (bestCounter.getCount() - secondCount >= this.minScore)
                retVal.setCount(bestCounter.getKey(), this.getWeight(len, covg), bestCounter.getRoleScores());
        }
        return retVal;
    }

}
