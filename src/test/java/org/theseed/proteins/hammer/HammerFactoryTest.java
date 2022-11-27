/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.DnaKmers;

/**
 * @author Bruce Parrello
 *
 */
class HammerFactoryTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HammerFactoryTest.class);

    @Test
    void testHammerFactory() throws IOException {
        File gFile = new File("data", "1002870.3.gto");
        DnaKmers.setKmerSize(20);
        Genome genome = new Genome(gFile);
        File rFile = new File("data", "roles.for.hammers");
        RoleMap roleMap = RoleMap.load(rFile);
        GenomeHammerFactory factory = new GenomeHammerFactory(genome, roleMap);
        // Verify the hammers.
        for (Feature feat : genome.getPegs()) {
            if (feat.isInteresting(roleMap)) {
                final String fid = feat.getId();
                log.info("Checking kmers for {}.", fid);
                var kmerStream = new SequenceKmerIterable(genome.getDna(fid), 20);
                for (String kmer : kmerStream)
                    assertThat(fid, factory.isPresent(kmer), equalTo(true));
            }
        }
    }

}
