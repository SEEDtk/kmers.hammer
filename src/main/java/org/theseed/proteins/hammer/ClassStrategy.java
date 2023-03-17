/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Collection;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;

/**
 * This object implements a strategy for classifying samples using hammers.  The main engine takes as input
 * a sequence and a list of hits, then returns a weight map for the sequence.  The weight map could indicate
 * a single group, multiple groups, or no group at all, but each group is associated with a score.  The
 * resulting weight map could be used to choose a bin or to accumulate scores for the whole sample.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ClassStrategy {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ClassStrategy.class);
    /** TRUE if the counts should be scaled, else FALSE */
    private boolean scaled;
    /** scale factor for weighting, based on the typical read length */
    protected static final double SCALE_FACTOR = 180.0;

    /**
     * This interface defines the methods that must be supported by any command processor that uses this object,
     * so that the subclasses can access the necessary tuning parameters.
     */
    public interface IParms {

        /**
         * @return the minimum threshold score required to be considered a dominant group for a sequence
         */
        public double getMinDiff();

        /**
         * @return TRUE to scale the counts by length and coverage, else FALSE
         */
        public boolean isScaled();

        /**
         * @return the counting method
         */
        public HammerDb.Method getMethod();

    }

    /**
     * This enumerates the types of strategies supported by the subclasses.
     */
    public static enum Type {
        /** report the raw number of hits on each sequence for all groups */
        HITS {
            @Override
            public ClassStrategy create(IParms processor) {
                return new HitsClassStrategy(processor);
            }
        },
        /** choose a single group to represent each sequence */
        REGIONS {
            @Override
            public ClassStrategy create(IParms processor) {
                return new RegionsClassStrategy(processor);
            }
        };

        /**
         * Create a strategy object of this type.
         *
         * @param processor		controlling command processor
         */
        public abstract ClassStrategy create(IParms processor);

    }

    /**
     * Construct a new strategy object for the specified command processor.
     *
     * @param processor		controlling command processor
     */
    public ClassStrategy(IParms processor) {
        this.scaled = processor.isScaled();
    }

    /**
     * Score a single sequence and return its weight map.
     *
     * @param hits		set of hits for the sequence
     * @param len		length of the sequence
     * @param covg		coverage of the sequence
     *
     * @return a weight map containing the score for each represented group in the sequence
     */
    public abstract WeightMap computeScores(Collection<HammerDb.Hit> hits, int len, double covg);

    /**
     * @return the weight to use for a specified length and coverage
     *
     * @param len	length of the sequence
     * @param covg	coverage of the sequence
     */
    public double getWeight(int len, double covg) {
        return (this.scaled ? len * covg / SCALE_FACTOR : 1.0);
    }

}
