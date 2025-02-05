package org.theseed.genome;

import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.proteins.DnaTranslator;
import org.theseed.proteins.RoleMap;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

public class GenomeUniSeqTest {

    /** kmer size for testing */
    private static final int KMER_SIZE = 20;
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeUniSeqTest.class);

    @Test
    public void testBadChars() throws IOException {
        GenomeUniSequences.setKmerSize(KMER_SIZE);
        File gtoFile = new File("data/seq_test", "28181.21.gto");
        RoleMap roles = RoleMap.load(new File("data", "roles.for.hammers"));
        Genome genome = new Genome(gtoFile);
        GenomeUniSequences seqs = new GenomeUniSequences(genome, roles);
        var kmerMap = seqs.getKmerMap();
        for (String kmer : kmerMap.keySet()) {
            int bad = StringUtils.indexOfAnyBut(kmer, "acgt");
            assertThat(kmer, bad, lessThan(0));
        }
    }

    @Test
    public void testSimple() throws IOException {
        GenomeUniSequences.setKmerSize(KMER_SIZE);
        File gtoFile = new File("data/seq_test", "171101.6.gto");
        RoleMap roles = RoleMap.load(new File("data", "roles.for.hammers"));
        assertThat(roles.size(), equalTo(5));
        Genome genome = new Genome(gtoFile);
        GenomeUniSequences seqs = new GenomeUniSequences(genome, roles);
        assertThat(seqs.getGenomeId(), equalTo("171101.6"));
        // Verify that we have the correct sequences.
        var seqMap = seqs.getSequenceMap();
        assertThat(seqMap.keySet(), containsInAnyOrder("fig|171101.6.peg.1917", "fig|171101.6.peg.1194",
                "fig|171101.6.peg.568", "fig|171101.6.peg.1064", "fig|171101.6.peg.559"));
        // Verify that we have the correct DNA.
        DnaTranslator xlate = new DnaTranslator(11);
        for (Map.Entry<String, String> seqEntry : seqMap.entrySet()) {
            String fid = seqEntry.getKey();
            Feature feat = genome.getFeature(fid);
            String prot = feat.getProteinTranslation();
            String dna = seqEntry.getValue();
            String prot2 = xlate.pegTranslate(dna);
            assertThat(fid, prot.toUpperCase() + "*", equalTo(prot2.toUpperCase()));
        }
        // Verify that the kmer map has the correct kmers.
        Map<String, String> kmerMap = seqs.getKmerMap();
        log.info("{} kmers found in {}.", kmerMap.size(), seqs);
        for (Map.Entry<String, String> seqEntry : seqMap.entrySet()) {
            String seq = seqEntry.getValue();
            String fid = seqEntry.getKey();
            final int n = seq.length() - KMER_SIZE;
            for (int i = 0; i <= n; i++) {
                String kmer = seq.substring(i, i+KMER_SIZE);
                assertThat(kmer, kmerMap.get(kmer), equalTo(fid));
            }
        }
        // Get a close genome.
        File otherFile = new File("data/seq_test", "1313.9775.gto");
        var genome2 = new Genome(otherFile);
        var seqs2 = new GenomeUniSequences(genome2, roles);
        seqs2.processKmers(kmerMap);
        Set<String> kmersLeft = kmerMap.keySet();
        log.info("{} kmers left after reduction.", kmersLeft.size());
        // Verify we removed the kmers.
        for (String seq : seqs2.getSequenceMap().values()) {
            final int n = seq.length() - KMER_SIZE;
            for (int i = 0; i <= n; i++) {
                var kmer1 = seq.substring(i, i+KMER_SIZE);
                assertThat(kmer1, not(in(kmersLeft)));
                assertThat(Contig.reverse(kmer1), not(in(kmersLeft)));
            }
        }
    }

}
