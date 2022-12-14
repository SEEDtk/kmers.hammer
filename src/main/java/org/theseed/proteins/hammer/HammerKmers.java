/**
 *
 */
package org.theseed.proteins.hammer;

import java.util.Set;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.sequence.DnaKmers;

/**
 * This is a special kmer object designed to manage hammer kmers.  When we are searching for hammers, we need
 * to remember the source feature ID, and we need the forward and reverse kmers stored separately.  We also
 * don't worry about distance or similarity or anything like that.  All we want to know is whether an incoming
 * kmer is found in this sequence.
 *
 * To save memory, the kmers are not stored in a set.  We scan every time.
 *
 * @author Bruce Parrello
 *
 */
public class HammerKmers {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerKmers.class);

    /** feature ID for the originating sequence */
    private String fid;
    /** sequence in the forward direction */
    private String forward;
    /** sequence on the reverse strand */
    private String reverse;
    /** kmer size to use */
    private static int K = 20;

    /**
     * This is a simple enum that defines the mode for a hammer-kmer check.
     */
    public static enum Mode {
        /** worthiness mode:  check forward only */
        WORTHINESS {
            @Override
            public boolean check(String hammer, HammerKmers kmerSet) {
                return kmerSet.isInSequence(hammer);
            }
        },
        /** precision mode:  check both directions */
        PRECISION {
            @Override
            public boolean check(String hammer, HammerKmers kmerSet) {
                return kmerSet.isInCommon(hammer);
            }
        };

        /**
         * Check a kmer set for a match in the current mode.
         *
         * @param hammer	hammer for which to look
         * @param kmerSet	kmer set to check
         *
         * @return TRUE if the hammer was found, else FALSE
         */
        public abstract boolean check(String hammer, HammerKmers kmerSet);

    }

    /**
     * Specify the global kmer size.
     *
     * @param kSize		kmer size to use for hammers
     *
     */
    public static void setKmerSize(int kSize) {
        K = kSize;
    }

    /**
     * @return the global kmer size
     */
    public static int getKmerSize() {
        return K;
    }

    /**
     * Construct a hammer-kmer object for a DNA sequence.
     *
     * @param seqId		ID of the sequence
     * @param seq		DNA string for the sequence
     */
    public HammerKmers(String seqId, String seq) {
        this.fid = seqId;
        // Store the forward sequence.
        this.forward = seq.toLowerCase();
        // Store the reverse sequence.
        this.reverse = Contig.reverse(this.forward);
    }

    /**
     * Return a set of all the kmers in the specified DNA sequence.
     *
     * @param seq		DNA sequence to process for kmers
     *
     * @return the set of kmers found
     */
    private static Set<String> collectKmers(String seq) {
        var iter = new SequenceKmerIterable(seq, K);
        Set<String> retVal = iter.stream().filter(x -> DnaKmers.isClean(x)).collect(Collectors.toSet());
        return retVal;
    }

    /**
     * @return the set of forward kmers
     */
    public Set<String> getKmers() {
        return collectKmers(this.forward);
    }

    /**
     * @return the feature ID for this kmer set
     */
    public String getFid() {
        return this.fid;
    }

    /**
     * @return TRUE if a kmer is found in the forward set
     */
    public boolean isInSequence(String kmer) {
        return this.forward.contains(kmer);
    }

    /**
     * @return TRUE if a kmer is found in either set
     */
    public boolean isInCommon(String kmer) {
        return this.forward.contains(kmer) || this.reverse.contains(kmer);
    }

    @Override
    public int hashCode() {
        int result = ((this.fid == null) ? 0 : this.fid.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof HammerKmers)) {
            return false;
        }
        HammerKmers other = (HammerKmers) obj;
        if (this.fid == null) {
            if (other.fid != null) {
                return false;
            }
        } else if (!this.fid.equals(other.fid)) {
            return false;
        }
        return true;
    }

}
