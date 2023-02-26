/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.Set;

import org.junit.jupiter.api.Test;
import org.theseed.genome.Contig;
import org.theseed.sequence.DnaKmers;

/**
 * @author Bruce Parrello
 *
 */
class TestHammerKmers {

    private static final String seq1 =
            "gtgacagatgtgcaagtgtccaatggctatgtgttgcacattggtttcctcaagtacggtacayytacgtgtggatgaccaggtgartggtgaactacgacgaggcgncgtcgtcgtccgctg";


    @Test
    void testHammerKmers() {
        HammerKmers testH = new HammerKmers("fig|10.peg.20", seq1);
        assertThat(testH.getFid(), equalTo("fig|10.peg.20"));
        assertThat(testH.isInSequence("gtccaatggctatgtgttgc"), equalTo(true));
        assertThat(testH.isInCommon("gcaacacatagccattggac"), equalTo(true));
        assertThat(testH.isInSequence("gcaacacatagccattggac"), equalTo(false));
        // Acid test.
        final int k = HammerKmers.getKmerSize();
        final int n = seq1.length() - k;
        for (int i = 0; i < n; i++) {
            String kmer = seq1.substring(i, i + k);
            if (DnaKmers.isClean(kmer)) {
                assertThat(kmer, testH.isInSequence(kmer), equalTo(true));
                assertThat(kmer, testH.isInCommon(kmer), equalTo(true));
                String rev = Contig.reverse(kmer);
                assertThat(rev, testH.isInCommon(rev), equalTo(true));
            }
        }
        Set<String> kmers = testH.getKmers();
        for (String kmer : kmers) {
            assertThat(kmer, DnaKmers.isClean(kmer), equalTo(true));
        }
    }

}
