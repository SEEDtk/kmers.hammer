/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 *
 * This is a slight variant of the score map that only stores the maximum role count.
 * *
 * @author Bruce Parrello
 *
 */
public class SummaryMap {

    /** underlying hash map */
    private HashMap<String, Count> map;
    /** sorter for key-sorting */
    private static final Comparator<Count> KEY_SORTER = new KeySorter();

    /**
     * Utility class to encapsulate the count.
     */
    public class Count implements Comparable<Count> {

        /** representative genome ID */
        private String key;
        /** number of roles found */
        private int roleCount;
        /** current count value */
        private double num;

        /**
         * Create a new count.
         *
         * @param key	key being counted
         */
        private Count(String key) {
            this.key = key;
            this.num = 0.0;
            this.roleCount = 0;
        }

        /**
         * @return the key counted by this counter.
         */
        public String getKey() {
            return this.key;
        }

        /**
         * @return the value counted for this key
         */
        public double getCount() {
            return this.num;
        }

        /**
         * @return the number of roles found for this key
         */
        public int getRoleCount() {
            return this.roleCount;
        }

        /**
         * Merge another count into this count.
         *
         * @param other		other counter to add
         */
        protected void merge(SummaryMap.Count other) {
            this.num += other.num;
            if (other.roleCount > this.roleCount)
                this.roleCount = other.roleCount;
        }

        /**
         * The sort order is highest count to lowest count, then from most roles to fewest roles.
         * If counts and role-set sizes are equal but the keys are not, we need to compare different.
         * In that case, we sort by the string representation.
         */
        @Override
        public int compareTo(SummaryMap.Count o) {
            // Note we compare the counts in reverse.
            int retVal = Double.compare(o.num, this.num);
            if (retVal == 0) {
                retVal = o.getRoleCount() - this.getRoleCount();
                if (retVal == 0)
                    retVal = this.key.toString().compareTo(o.key.toString());
            }
            return retVal;
        }

        /**
         * The string representation includes the key and the count.
         */
        @Override
        public String toString() {
            String retVal = String.format("%s (weight %s)", this.key.toString(),
                        this.num);
            return retVal;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + getEnclosingInstance().hashCode();
            result = prime * result + ((key == null) ? 0 : key.hashCode());
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            Count other = (Count) obj;
            if (!getEnclosingInstance().equals(other.getEnclosingInstance()))
                return false;
            if (key == null) {
                if (other.key != null)
                    return false;
            } else if (!key.equals(other.key))
                return false;
            if (num != other.num)
                return false;
            return true;
        }

        private SummaryMap getEnclosingInstance() {
            return SummaryMap.this;
        }

    }

    /**
     * This is a comparator that sorts the counts by key.
     */
    public static class KeySorter implements Comparator<Count> {

        @Override
        public int compare(Count o1, Count o2) {
            return o1.key.compareTo(o2.key);
        }

    }


    /**
     * Create a blank weighted-counting map.
     */
    public SummaryMap() {
        this.map = new HashMap<String, Count>();
    }

    /*
     * Create an empty weighted-counting map with the specified initial capacity.
     *
     * @param capacity	number of keys expected
     */
    public SummaryMap(int capacity) {
        this.map = new HashMap<String, Count>(capacity * 4 / 3 + 1);
    }

    /**
     * @return the weighted count for a given key (which is 0 if the key is unknown)
     *
     * @param key	key for the counter of interest
     */
    public double getCount(String key) {
        double retVal = 0;
        Count found = this.map.get(key);
        if (found != null) {
            retVal = found.num;
        }
        return retVal;
    }

    /**
     * Erase all the counts in this map without deleting any keys.
     */
    public void clear() {
        for (Count count : this.map.values()) {
            count.num = 0.0;
            count.roleCount = 0;
        }
    }

    /**
     * @return the counter object for a specified key.  If none exists, one will be created.
     *
     * @param key	key for the counter of interest
     */
    protected Count getCounter(String key) {
        Count retVal = this.map.computeIfAbsent(key, x -> new Count(x));
        return retVal;
    }

    /**
     * @return the counter object for the specified key, or NULL if the key has not been counted
     *
     * @param key	key of interest
     */
    public Count findCounter(String key) {
        return this.map.get(key);
    }

    /**
     * @return	a sorted collection of all the keys in this object
     */
    public SortedSet<String> keys() {
        TreeSet<String> retVal = new TreeSet<String>(this.map.keySet());
        return retVal;
    }

    /**
     * @return the number of keys in this map
     */
    public int size() {
        return this.map.size();
    }

    /**
     * @return a list of all the counts in this object, sorted from highest to lowest
     */
    public List<Count> sortedCounts() {
        ArrayList<Count> retVal = new ArrayList<Count>(this.map.values());
        retVal.sort(null);
        return retVal;
    }

    /**
     * @return a list of all the counts in the object, sorted by key
     */
    public List<Count> keyedCounts() {
        ArrayList<Count> retVal = new ArrayList<Count>(this.map.values());
        retVal.sort(KEY_SORTER);
        return retVal;
    }

    /**
     * @return an unordered list of all the counts
     */
    public Collection<Count> counts() {
        return this.map.values();
    }

    /**
     * Erase all classes and counts from this map.
     */
    public void deleteAll() {
        this.map.clear();
    }

    /**
     * Set the count to a specific value.
     *
     * @param key		key whose count is to be set
     * @param newValue	value of the new weighted count
     * @param roleCount number of roles to put into the count
     */
    public void setCount(String key, double newValue, int roleCount) {
       Count myCount = this.getCounter(key);
       myCount.num = newValue;
       myCount.roleCount = roleCount;
    }

    /**
     * @return the sum of all the counts in this map
     */
    public double sum() {
        return this.map.values().stream().mapToDouble(x -> x.getCount()).sum();
    }

    /**
     * Remove a count from the map.
     *
     * @param key	key of the count
     */
    public void remove(String key) {
        this.map.remove(key);
    }

    /**
     * @return the counter entry with the highest count, or NULL if there are no counts
]	 */
    public Count getBestEntry() {
        Count retVal = null;
        // We get the minimum, since counts sort with the highest first.
        if (this.map.size() > 0)
            retVal = Collections.min(this.map.values());
        return retVal;
    }

    /**
     * Add a new count value for the specified key.
     *
     * @param key			key for the count
     * @param count			value to add
     * @param roleCount		number of roles
     */
    public void count(String key, double count, int roleCount) {
        Count myCount = this.getCounter(key);
        myCount.num += count;
        if (roleCount > myCount.roleCount)
            myCount.roleCount = roleCount;
    }

}
