/**
 *
 */
package org.theseed.binning;

import java.io.IOException;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.counters.CountMap;

/**
 * This is the base class for deciding which bin should contain each contig based on the hammer hits.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerBinningRule {

    /**
     * This interface must be supported by any command processor used to construct a binning rule.
     */
    public interface IParms {

        /**
         * @return the minimum difference between the highest value and the second-highest value
         */
        public int getMinDiff();

    }

    /**
     * This enum defines the different binning rule types.
     */
    public static enum Type {
        /** choose the bin with the most hits */
        MAX {
            @Override
            public HammerBinningRule create(IParms processor) throws ParseFailureException {
                return new MaxHammerBinningRule(processor);
            }
        };

        public abstract HammerBinningRule create(IParms processor) throws IOException, ParseFailureException;
    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerBinningRule.class);

    /**
     * Choose the bin for a contig.
     *
     * @param counts	count map of the number of hits per bin
     *
     * @return the bin ID, or NULL if the contig should be skipped
     */
    public abstract String choose(CountMap<String> counts);


}
