/**
 *
 */
package org.theseed.sequence.seeds.filters;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.sequence.Sequence;

/**
 * This filter is used to determine which features should be included when processing a finder
 * role file to produce hammers.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerFeatureFilter {

     /**
     * This interface provides access to the parameters that must be provided by
     * a command processor that supports feature filters.
     */
    public interface IParms {

    }

    /**
     * This enum describes the different types of feature filters.
     */
    public static enum Type {
        /** keep all features */
        ALL {
            @Override
            public HammerFeatureFilter create(IParms processor) {
                return new HammerFeatureFilter.All(processor);
            }
        },
        /** only keep features with single roles */
        SIMPLE {
            @Override
            public HammerFeatureFilter create(IParms processor) {
                return new HammerFeatureFilter.Simple(processor);
            }
        };

        /**
         * @return a hammer-feature filter for use by the specified command processor
         *
         * @param processor		controlling command processor
         */
        public abstract HammerFeatureFilter create(IParms processor);

    }

    /**
     * Construct a hammer feature filter for the specified command processor.
     *
     * @param processor		controlling command processor
     */
    public HammerFeatureFilter(IParms processor) {
    }

    /**
     * @return TRUE if the specified sequence should be used during hammer generation, else FALSE
     *
     * @param seq	sequence (id, comment, DNA) to check
     */
    public abstract boolean check(Sequence seq);

    /**
      * This filter does not filter out anything.
      */
     public static class All extends HammerFeatureFilter {

         public All(IParms processor) {
             super(processor);
         }

         @Override
         public boolean check(Sequence seq) {
             return true;
         }

     }

     /**
      * This filter removes features with compound roles.
      */
     public static class Simple extends HammerFeatureFilter {

         public Simple(IParms processor) {
             super(processor);
         }

         @Override
         public boolean check(Sequence seq) {
             // Get the comment and split it into fields.
             String[] parts = StringUtils.split(seq.getComment(), '\t');
             String roles[] = Feature.rolesOfFunction(parts[2]);
             return (roles.length == 1);
         }

     }




}
