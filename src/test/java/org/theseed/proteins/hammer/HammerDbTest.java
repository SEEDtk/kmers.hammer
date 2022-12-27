/**
 *
 */
package org.theseed.proteins.hammer;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.theseed.counters.CountMap;
import org.theseed.genome.Genome;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.sqlite.SqliteDbConnection;
import org.theseed.locations.Location;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ParseFailureException;

/**
 * @author Bruce Parrello
 *
 */
class HammerDbTest {

    private static final File LOAD_FILE = new File("data", "hammers200.tbl");

    @Test
    void testHashHammerDb() throws IOException, ParseFailureException {
        HashHammerDb hammers = new HashHammerDb(LOAD_FILE);
        hammerDbTester(hammers);
    }

    @Test
    void testSqlHammerDb() throws IOException, ParseFailureException, SQLException {
        File dbFile = new File("data", "hammerdb.ser");
        if (dbFile.exists())
            FileUtils.forceDelete(dbFile);
        DbConnection db = new SqliteDbConnection(dbFile);
        File initFile = new File("data", "hammerdb.sql");
        db.scriptUpdate(initFile);
        SqlHammerDb hammers = new SqlHammerDb(db, 500);
        hammers.load(LOAD_FILE);
        // Now the database is loaded, so we run the tests.
        hammerDbTester(hammers);
    }

    /**
     * This tests a hammer database.  The database must be loaded before calling from
     * hammers200.tbl.
     *
     * @param hammers	loaded hammer database
     *
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected void hammerDbTester(HammerDb hammers) throws FileNotFoundException, IOException {
        List<Sequence> seqs = new ArrayList<Sequence>();
        try (FastaInputStream seqsIn = new FastaInputStream(new File("data", "hammer.test.fa"))) {
            for (Sequence seq : seqsIn)
                seqs.add(seq);
        }
        CountMap<String> counts = hammers.findClosest(seqs);
        assertThat(counts.size(), equalTo(3));
        List<CountMap<String>.Count> results = counts.sortedCounts();
        assertThat(results.get(0).getKey(), equalTo("565575.4"));
        assertThat(results.get(1).getKey(), equalTo("1397.4"));
        assertThat(results.get(2).getKey(), equalTo("1278308.3"));
        File gFile = new File("data", "2485170.3.gto");
        List<Sequence> gSeqs = new Genome(gFile).getSequences();
        seqs.addAll(gSeqs);
        counts = hammers.findClosest(seqs);
        results = counts.sortedCounts();
        assertThat(results.size(), equalTo(4));
        assertThat(results.get(0).getKey(), equalTo("2485170.3"));
        File hitFile = new File("data", "hammer.hit.test.fa");
        seqs = FastaInputStream.readAll(hitFile);
        var hits = hammers.findHits(seqs);
        assertThat(hits.size(), equalTo(7));
        Set<String> fidsHit = hits.stream().map(x -> x.getFid()).collect(Collectors.toSet());
        assertThat(fidsHit, containsInAnyOrder("fig|1278308.3.peg.2084", "fig|1278308.3.peg.2569", "fig|1397.4.peg.5364"));
        for (var hit : hits) {
            String fid = hit.getFid();
            switch (fid) {
            case "fig|1278308.3.peg.2084" :
                assertThat(fid, hit.getLoc(), anyOf(equalTo(Location.create("seq2", 30, 11)),
                        equalTo(Location.create("seq2", 31, 12))));
                assertThat(fid, hit.getStrength(), closeTo(0.6, 0.05));
                break;
            case "fig|1278308.3.peg.2569" :
                assertThat(fid, hit.getLoc(), anyOf(equalTo(Location.create("seq2", 49, 68)),
                        equalTo(Location.create("seq1", 129, 110)),
                        equalTo(Location.create("seq2", 50, 69)),
                        equalTo(Location.create("seq2", 51, 70))));
                assertThat(fid, hit.getStrength(), closeTo(0.7, 0.05));
                break;
            case "fig|1397.4.peg.5364" :
                assertThat(fid, hit.getLoc(), equalTo(Location.create("seq4", 42, 61)));
                assertThat(fid, hit.getStrength(), closeTo(0.8, 0.05));
                break;
            }
        }

    }

}
