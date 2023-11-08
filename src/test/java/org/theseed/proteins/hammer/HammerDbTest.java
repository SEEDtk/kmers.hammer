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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.sqlite.SqliteDbConnection;
import org.theseed.locations.Location;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.utils.ParseFailureException;

/**
 * @author Bruce Parrello
 *
 */
class HammerDbTest {

    private static final File LOAD_FILE = new File("data", "hammers200.tbl");
    protected static Logger log = LoggerFactory.getLogger(HammerDbTest.class);


    @Test
    void testHashHammerDb() throws IOException, ParseFailureException {
        HashHammerDb hammers = new HashHammerDb(LOAD_FILE);
        hammerDbTester(hammers);
        assertThat(hammers.getLoadFile(), equalTo(LOAD_FILE));
    }

    @Test
    void testQualityStuff() throws IOException, ParseFailureException {
        double q = SeqRead.qualChance("FFFFFFFFFFFFFFFFFFFF", 0, 20);
        log.info("Normal quality value is {}.", q);
        HashHammerDb hammers = new HashHammerDb(LOAD_FILE);
        SeqRead.Part part1 = new SeqRead.Part("part1", "xxxxxxxgaggtcgacaacgacatcgcxxxxxxxxgcgatgtcgttgtcgacctcxxxxxxxx",
                                                       "FFFFFFF,,,,,,,,,,,,,,,,,,,,FFFFFFFF::::::::::::::::::::FFFFFFFF");
        SeqRead.Part part2 = new SeqRead.Part("part2", "xxxxxxxatcgcggcggctgccgacgcxxxxxxxxgcgatgtcgttgtcgacctcxxxxxxxx",
                                                          "FFFFFFF,,,,,,,,,,,,,,,,,,,,FFFFFFFF::::::::::::::::::::FFFFFFFF");
        List<SeqRead.Part> parts = List.of(part1, part2);
        SortedSet<HammerDb.Hit> hits = hammers.findHits(parts, 0.0);
        Iterator<HammerDb.Hit> iter = hits.iterator();
        var hit1 = iter.next();
        Location loc = hit1.getLoc();
        assertThat(loc.getContigId(), equalTo("part1"));
        assertThat(loc.getDna(part1.getSequence()), equalTo("gaggtcgacaacgacatcgc"));
        assertThat(part1.getQual().substring(loc.getLeft() - 1, loc.getRight()), equalTo(",,,,,,,,,,,,,,,,,,,,"));
        hit1 = iter.next();
        loc = hit1.getLoc();
        assertThat(loc.getContigId(), equalTo("part1"));
        assertThat(loc.getDna(part1.getSequence()), equalTo("gaggtcgacaacgacatcgc"));
        assertThat(part1.getQual().substring(loc.getLeft() - 1, loc.getRight()), equalTo("::::::::::::::::::::"));
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
     * This tests a hammer database.  The database must be loaded from hammers200.tbl before calling.
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
        WeightMap counts = hammers.findClosest(seqs);
        assertThat(counts.size(), equalTo(3));
        List<WeightMap.Count> results = counts.sortedCounts();
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
        // Test the hammers-in-genome method.  We need to re-read the hammer file to get a list.
        Set<String> found = new HashSet<String>(8800);
        try (TabbedLineReader dbStream = new TabbedLineReader(LOAD_FILE)) {
            for (var line : dbStream) {
                String fid = line.get(1);
                if (fid.startsWith("fig|1278308.3.peg"))
                    found.add(line.get(0));
            }
        }
        // Now compare the lists.
        var gHammers = hammers.findGenomeHammers("1278308.3");
        assertThat(gHammers.size(), equalTo(found.size()));
        for (String hammer : gHammers.keySet())
            assertThat(hammer, in(found));
    }



}
