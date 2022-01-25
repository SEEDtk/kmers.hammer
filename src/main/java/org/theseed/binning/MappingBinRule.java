/**
 *
 */
package org.theseed.binning;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import org.theseed.sequence.Sequence;

/**
 * This object uses a hash of contig IDs to bin IDs in order to determine which bin contains
 * a contig.
 *
 * @author Bruce Parrello
 *
 */
public class MappingBinRule implements BinBuilder.IRule {

    // FIELDS
    /** map of contig IDs to bin IDs */
    private Map<String, String> contigMap;
    /** minimum contig coverage */
    private double minCovg;
    /** minimum contig length */
    private int minLen;
    /** default minimum length */
    public static final int DEFAULT_MIN_LEN = 400;
    /** default minimum coverage */
    public static final double DEFAULT_MIN_COVERAGE = 4.0;

    /**
     * Construct a map-based binning rule.
     *
     * @param map	map of contig IDs to bin IDs
     */
    public MappingBinRule(Map<String, String> map) {
        this.contigMap = map;
        this.minCovg = DEFAULT_MIN_COVERAGE;
        this.minLen = DEFAULT_MIN_LEN;
    }

    @Override
    public Set<String> getBinIds() {
        Set<String> retVal = new HashSet<String>(this.contigMap.values());
        return retVal;
    }

    @Override
    public String choose(Sequence seq) {
        // Filter for the contig coverage and length.
        String retVal = null;
        if (BinBuilder.contigFilter(seq, this.minCovg, this.minLen))
            retVal = this.contigMap.get(seq.getLabel());
        return retVal;
    }

    /**
     * Specify the minimum acceptable contig coverage.
     *
     * @param minCovg 	the minCovg to set
     */
    public void setMinCovg(double minCovg) {
        this.minCovg = minCovg;
    }

    /**
     * Specify the minimum acceptable contig length.
     *
     * @param minLen 	the minLen to set
     */
    public void setMinLen(int minLen) {
        this.minLen = minLen;
    }

}
