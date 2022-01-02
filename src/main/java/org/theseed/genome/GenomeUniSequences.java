/**
 *
 */
package org.theseed.genome;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.theseed.locations.Location;
import org.theseed.proteins.RoleMap;

/**
 * This object tracks the roles of interest for a genome.  It contains the genome ID plus
 * a map from feature IDs to universal role DNA sequences.  It is the primary workhorse of
 * the hammer search.
 *
 * Because we expect the number of sequences to be small, we use a tree map to conserve memory.
 * The objects are compared by genome ID.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeUniSequences implements Comparable<GenomeUniSequences> {

    // FIELDS
    /** ID of the genome of interest */
    private String genomeId;
    /** map of feature IDs to role sequences */
    private Map<String, String> sequenceMap;
    /** DNA kmer size */
    private static int kmerSize = 20;
    /** total base pairs in all sequences */
    private int seqLen;

    /**
     * Construct a role sequences object from a role map and a genome.
     *
     * @param genome	source genome
     * @param roleMap	role definition map for the universal roles to use
     */
    public GenomeUniSequences(Genome genome, RoleMap roleMap) {
        // Save the genome ID.
        this.genomeId = genome.getId();
        // Create the sequence map.
        this.sequenceMap = new TreeMap<String, String>();
        // Denote no sequences were found.
        this.seqLen = 0;
        // Loop through the pegs, searching for roles of interest.
        for (Feature feat : genome.getPegs()) {
            if (feat.getUsefulRoles(roleMap).size() > 0) {
                Location loc = feat.getLocation();
                this.sequenceMap.put(feat.getId(), genome.getDna(loc));
                this.seqLen += loc.getLength();
            }
        }
    }

    /**
     * Specify a new kmer size.
     *
     * @param k		new kmer size
     */
    public static void setKmerSize(int k) {
        kmerSize = k;
    }

    /**
     * Compute the kmer map for this genome.
     *
     * @return a map from forward DNA kmers to feature IDs
     */
    public Map<String, String> getKmerMap() {
        // Compute the size for the hash map.
        int hashLen = this.seqLen * 4 / 3;
        // Create the map.
        var retVal = new HashMap<String, String>(hashLen);
        for (Map.Entry<String, String> seqEntry : this.sequenceMap.entrySet()) {
            String fid = seqEntry.getKey();
            List<String> parts = Contig.cleanParts(seqEntry.getValue());
            for (String seq : parts) {
                final int n = seq.length() - kmerSize;
                for (int i = 0; i <= n; i++)
                    retVal.put(seq.substring(i, i+kmerSize), fid);
            }
        }
        return retVal;
    }

    /**
     * Eliminate all kmers found in this genome's sequences from a kmer map.
     *
     * @param kmerMap	map of kmers to feature IDs (to be modified)
     */
    public void processKmers(Map<String, String> kmerMap) {
        // Loop through the sequences.
        for (String seq : this.sequenceMap.values()) {
            // Eliminate all the forward kmers.
            this.processSeqKmers(kmerMap, seq);
            // Reverse-complement and eliminate all the backward kmers.
            String revComp = Contig.reverse(seq);
            this.processSeqKmers(kmerMap, revComp);
        }
    }

    /**
     * Eliminate all the kmers in the specified sequence from the incoming
     * kmer map.
     *
     * @param kmerMap	map of kmers to feature IDs (to be modified)
     * @param seq		DNA sequence to scan
     */
    private void processSeqKmers(Map<String, String> kmerMap, String seq) {
        int n = seq.length() - kmerSize;
        for (int i = 0; i <= n; i++)
            kmerMap.remove(seq.substring(i, i + kmerSize));
    }

    /**
     * @return the sequence map for this genome
     */
    protected Map<String, String> getSequenceMap() {
        return this.sequenceMap;
    }

    /**
     * @return the ID of this genome
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    @Override
    public int compareTo(GenomeUniSequences o) {
        if (o == null)
            return 1;
        else
            return this.genomeId.compareTo(o.genomeId);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((this.genomeId == null) ? 0 : this.genomeId.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof GenomeUniSequences)) {
            return false;
        }
        GenomeUniSequences other = (GenomeUniSequences) obj;
        if (this.genomeId == null) {
            if (other.genomeId != null) {
                return false;
            }
        } else if (!this.genomeId.equals(other.genomeId)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return String.format("%s (%d sequences, %d base pairs)",
                this.genomeId, this.sequenceMap.size(), this.seqLen);
    }

}
