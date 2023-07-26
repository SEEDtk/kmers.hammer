/**
 *
 */
package org.theseed.sequence.seeds.filters;

import java.util.Map;
import java.util.TreeMap;

/**
 * This class examines a list of duplicate features for a single role in a single genome and determines
 * which ones to keep.  It is used when computing hammers from a finder.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerDupFilter {

    // FIELDS
    /** empty map to use when returning nothing */
    private static final Map<String, String> EMPTY_MAP = new TreeMap<String, String>();

    /**
     * This interface provides access to the parameters that must be provided by
     * a command processor that supports feature filters.
     */
    public interface IParms {

    }

    /**
     * This enum specifies the types of duplicate-feature filters supported.
     */
    public static enum Type {
        /** keep all duplicates */
        KEEP {
            @Override
            public HammerDupFilter create(IParms processor) {
                return new HammerDupFilter.Keep(processor);
            }
        },
        /** keep only the longest protein */
        LONGEST {
            @Override
            public HammerDupFilter create(IParms processor) {
                return new HammerDupFilter.Longest(processor);
            }
        },
        /** reject all duplicates */
        NONE {
            @Override
            public HammerDupFilter create(IParms processor) {
                return new HammerDupFilter.None(processor);
            }
        };

        /**
         * @return a duplicate-feature filter for hammer generation by the specified command processor
         *
         * @param processor		controlling command processor
         */
        public abstract HammerDupFilter create(IParms hammerFinderProcessor);

    }

    /**
     * Construct a duplicate-feature filter for a specified command processor.
     *
     * @param processor		controlling command processor
     */
    public HammerDupFilter(IParms processor) {
    }

    /**
     * @return a reduced feature map with the appropriate duplicates filtered out
     *
     * @param fidMap	input feature map, with all examples of a genome's features for a single role,
     * 					mapping feature IDs to DNA sequences
     */
    public abstract Map<String, String> filter(Map<String, String> fidMap);

    /**
     * This filter keeps everything.
     */
    public static class Keep extends HammerDupFilter {

        public Keep(IParms processor) {
            super(processor);
        }

        @Override
        public Map<String, String> filter(Map<String, String> fidMap) {
            return fidMap;
        }

    }


    /**
     * This filter keeps nothing.
     */
    public static class None extends HammerDupFilter {

        public None(IParms processor) {
            super(processor);
        }

        @Override
        public Map<String, String> filter(Map<String, String> fidMap) {
            Map<String, String> retVal;
            if (fidMap.size() == 1)
                retVal = fidMap;
            else
                retVal = EMPTY_MAP;
            return retVal;
        }

    }


    /**
     * This filter keeps the longest protein.
     */
    public static class Longest extends HammerDupFilter {

        public Longest(IParms processor) {
            super(processor);
        }

        @Override
        public Map<String, String> filter(Map<String, String> fidMap) {
            Map<String, String> retVal;
            if (fidMap.size() <= 1)
                retVal = fidMap;
            else {
                retVal = new TreeMap<String, String>();
                // Find the longest feature.
                int len = 0;
                String fid = null;
                for (var fidEntry : fidMap.entrySet()) {
                    int newLen = fidEntry.getValue().length();
                    if (newLen > len) {
                        fid = fidEntry.getKey();
                        len = newLen;
                    }
                }
                // Save it if we found one.
                if (fid != null)
                    retVal.put(fid, fidMap.get(fid));
            }
            return retVal;
        }

    }



}
