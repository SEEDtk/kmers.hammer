/**
 *
 */
package org.theseed.binning;

import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;

/**
 * This binning rule selects the bin with the highest number of hits, provided it is
 * sufficiently higher than the next-highest hit count.
 *
 * @author Bruce Parrello
 *
 */
public class MaxHammerBinningRule extends HammerBinningRule {

    // FIELDS
    /** minimum hit-count distance */
    private int minHitDiff;

    /**
     * Construct a max-hits hammer binning rule.
     *
     * @param processor		controlling command processor
     *
     * @throws ParseFailureException
     */
    public MaxHammerBinningRule(IParms processor) throws ParseFailureException {
        this.minHitDiff = processor.getMinDiff();
        if (this.minHitDiff < 1)
            throw new ParseFailureException("Minimum hit difference must be positive.");
    }

    @Override
    public String choose(CountMap<String> counts) {
        String retVal = null;
        var counters = counts.sortedCounts();
        int nCounts = counters.size();
        if (nCounts > 0) {
            // Get the number of hits for the best bin.
            int best = counters.get(0).getCount();
            // Note that if only one bin is hit, the second-best is treated as 0 hits.
            int diff = (nCounts == 1 ? best : best - counters.get(1).getCount());
            if (diff >= this.minHitDiff)
                retVal = counters.get(0).getKey();
        }
        return retVal;
    }

}
