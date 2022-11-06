/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * This object takes as input a genome, a role map, and a FASTA file of contigs.  The contigs belong to genomes, and the
 * comment field should contain a genome ID.  The goal will be to generate hammers for the genomes in question.  A kmer
 * is considered a hammer if it is from the coding region of one of the roles in the map, and it is not in any of the
 * other genome contigs.
 *
 * The basic strategy is to create a hash of kmers to feature IDs from the genome. The FASTA file is then read one
 * contig at a time, and any kmers found in a contig that does NOT belong to the input genome are deleted.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeHammerFactory {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeHammerFactory.class);
    /** map of kmers to feature IDs; an ambiguous kmer is mapped to an empty string */
    private ConcurrentMap<String, String> kmerMap;
    /** ID of the genome containing the hammers */
    private String genomeId;
    /** TRUE for parallel processing, else FALSE */
    private static boolean paraMode = true;

    /**
     * Specify whether or not we should use parallel processing.
     *
     * @param paraMode 	TRUE to turn on parallel processing, FALSE to turn it off
     */
    public static void setParaMode(boolean paraMode) {
        GenomeHammerFactory.paraMode = paraMode;
    }

    /**
     * Create a new hammer factory.  The specified genome is scanned to create the initial kmer hash.
     *
     * @param genome	source genome
     * @param roles		definition map for roles of interest
     */
    public GenomeHammerFactory(Genome genome, RoleMap roles) {
        this.kmerMap = new ConcurrentHashMap<String, String>(roles.size() * 1400);
        this.genomeId = genome.getId();
        // Loop through the genome, checking the pegs.
        int usefulCount = 0;
        for (Feature feat : genome.getPegs()) {
            if (feat.getUsefulRoles(roles).size() > 0) {
                usefulCount++;
                // Here we have a feature of interest.  Extract its kmers in each direction.
                String fid = feat.getId();
                String fSeq = genome.getDna(fid);
                this.extractKmers(fSeq, fid);
                this.extractKmers(Contig.reverse(fSeq), fid);
            }
        }
        log.info("{} kmers found in {} useful features.", this.kmerMap.size(), usefulCount);
    }

    /**
     * Extract the kmers from the specified sequence (representing the feature with the given ID)
     * and update the kmer map.
     *
     * @param seq	sequence to scan
     * @param fid	ID of the feature containing the sequence
     */
    private void extractKmers(String seq, String fid) {
        var kmerStream = new SequenceKmerIterable(seq, DnaKmers.kmerSize());
        kmerStream.stream(paraMode).forEach(x -> this.kmerMap.put(x, fid));
    }

    /**
     * @return the feature ID for a kmer, or NULL if the kmer is not a hammer
     *
     * @param kmer	kmer to check
     */
    public String getFid(String kmer) {
        return this.kmerMap.get(kmer);
    }

    /**
     * @return TRUE if a kmer is present in the map, else FALSE
     *
     * @param kmer	kmer to check
     */
    public boolean isPresent(String kmer) {
        return this.kmerMap.containsKey(kmer);
    }

    /**
     * Read the FASTA file and delete non-discriminating kmers.
     *
     * @param file	FASTA file to read; the comments must contain the genome IDs
     *
     * @throws IOException
     */
    public void processFasta(File file) throws IOException {
        int oldCount = this.kmerMap.size();
        try (FastaInputStream fastaStream = new FastaInputStream(file)) {
            int count = 0;
            long timer = System.currentTimeMillis();
            for (Sequence seq : fastaStream) {
                if (! seq.getComment().contentEquals(this.genomeId)) {
                    var kmerStream = new SequenceKmerIterable(seq, DnaKmers.kmerSize());
                    kmerStream.stream(paraMode).forEach(x -> this.kmerMap.remove(x));
                    count++;
                    // Drop a progress message every 5 seconds or so.
                    if (log.isInfoEnabled() && System.currentTimeMillis() - timer >= 5000) {
                        timer = System.currentTimeMillis();
                        log.info("{} contigs processed.", count);
                    }
                }
            }
        }
        int newCount = this.kmerMap.size();
        log.info("{} kmers remaining after processing {}.  {} removed.", newCount, file, oldCount - newCount);
    }

    /**
     * Write the hammers to an output stream.
     *
     * @param writer	output stream for the hammers
     */
    public void dumpHammers(PrintWriter writer) {
        for (var kmerEntry : this.kmerMap.entrySet())
            writer.println(kmerEntry.getKey() + "\t" + kmerEntry.getValue());
    }

}
