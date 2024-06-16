/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.theseed.counters.CountMap;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.DbLoader;
import org.theseed.java.erdb.DbQuery;
import org.theseed.java.erdb.DbRecord;
import org.theseed.java.erdb.Relop;
import org.theseed.java.erdb.SqlBuffer;
import org.theseed.sequence.ISequence;
import org.theseed.sequence.KmerSeries;
import org.theseed.sequence.Sequence;
import org.theseed.stats.WeightMap;

/**
 * This is a hammer database stored in an SQL database.  It is slower than the in-memory
 * database, but is less stressful on the computer.
 *
 * @author Bruce Parrello
 *
 */
public class SqlHammerDb extends HammerDb {

    // FIELDS
    /** connection to the database */
    private DbConnection db;
    /** batch size to use */
    private int batchSize;
    /** size of a hash to hold a kmer batch */
    private int hashSize;

    /**
     * Construct an SQL hammer database for a specific command processor.
     *
     * @param processor		controlling command processor
     */
    public SqlHammerDb(HammerDb.IParms processor) {
        this.db = processor.getConnection();
        this.batchSize = processor.getBatchSize();
        this.init();
    }

    /**
     * Construct an SQL hammer database.
     *
     * @param db	database connection to use
     * @param b		batch size for queries
     */
    public SqlHammerDb(DbConnection db, int b) {
        this.db = db;
        this.batchSize = b;
        this.init();
    }

    /**
     * Initialize the object fields.
     */
    protected void init() {
        this.hashSize = this.batchSize * 4 / 3 + 1;
        // We need to compute the kmer size.  Get one kmer.
        SqlBuffer sizeQuery = new SqlBuffer(this.db).append("SELECT ").quote("Hammer", "hammer")
                .append(" FROM ").quote("Hammer").append(" LIMIT 1");
        try (PreparedStatement stmt = this.db.createStatement(sizeQuery)) {
            ResultSet results = stmt.executeQuery();
            if (! results.next())
                log.warn("No hammers found in hammer database.");
            else {
                String hammer = results.getString(1);
                final int kSize = hammer.length();
                this.setKmerSize(kSize);
                log.info("Kmer size is {}.", kSize);
            }
        } catch (SQLException e) {
            throw new RuntimeException("SQL error computing kmer size: " + e.toString());
        }
    }

    protected class Loader implements HammerDb.ILoader {

        /** database loader object for Hammer table */
        private DbLoader loader;
        /** controlling database transaction */
        private DbConnection.Transaction xact;

        /**
         * Construct the loader.
         */
        private Loader() {
            try {
                this.xact = SqlHammerDb.this.db.new Transaction();
                this.loader = DbLoader.batch(SqlHammerDb.this.db, "Hammer");
            } catch (SQLException e) {
                throw new RuntimeException("SQL Error: " + e.toString());
            }
        }

        @Override
        public void updateHammerMap(String fid, String roleId, String hammer, double str) {
            try {
                this.loader.set("hammer", hammer);
                this.loader.set("role", roleId);
                this.loader.set("fid", fid);
                this.loader.set("strength", str);
                this.loader.insert();
            } catch (SQLException e) {
                throw new RuntimeException("SQL Error: " + e.toString());
            }
        }

        @Override
        public void createEmptyMap(File inFile) throws SQLException {
            // Here we must erase the table if it has records in it.
            SqlBuffer buffer = SqlHammerDb.this.db.buildDeleteStmt("Hammer");
            try (PreparedStatement stmt = SqlHammerDb.this.db.createStatement(buffer)) {
                log.info("Clearing Hammer table prior to load.");
                stmt.execute();
            }
        }

        @Override
        public void close() {
            try {
                this.loader.close();
                this.xact.commit();
                this.xact.close();
            } catch (SQLException e) {
                throw new RuntimeException("SQL Error in cleanup: " + e.toString());
            }
        }

    }

    @Override
    protected ILoader getLoader() {
        return this.new Loader();
    }

    @Override
    protected void findClosestInternal(WeightMap map, Collection<Sequence> seqs, int kSize) {
        // For performance reason, we batch the kmer queries.  We will collect a batch of
        // kmers, ask for the feature IDs, and count the results.  This gets repeated
        // until we reach the end.
        try {
            // Get an iterator for the kmers.
            Iterable<String> kIter = KmerSeries.init(seqs, kSize);
            log.debug("Iterating through {} sequences.", seqs.size());
            // This will hold the current kmer batch.
            var kmers = new CountMap<String>();
            // Create a query.
            DbQuery query = null;
            try {
                 query = this.buildBatchQuery(this.batchSize);
                // Loop through the kmers, building a batch.
                for (String kmer : kIter) {
                    kmers.count(kmer);
                    if (kmers.size() == this.batchSize) {
                        // Now we run the query.
                        this.runCountQuery(query, kmers, map);
                        // Clear the kmer set for the next query.
                        kmers.deleteAll();
                    }
                }
            } finally {
                if (query != null)
                    query.close();
            }
            // Do we have a residual?
            if (kmers.size() > 0) {
                // Yes.  Create a mini-query for it.
                query = null;
                try {
                    query = this.buildBatchQuery(kmers.size());
                    this.runCountQuery(query, kmers, map);
                } finally {
                    if (query != null)
                        query.close();
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException("Find-closest query error: " + e.toString());
        }
    }

    /**
     * @return a batch hammer query for counts and hits
     *
     * @param size	number of kmers to be passed in
     *
     * @throws SQLException
     */
    protected DbQuery buildBatchQuery(int size) throws SQLException {
        DbQuery retVal = new DbQuery(this.db, "Hammer");
        retVal.select("Hammer", "fid", "hammer", "strength", "role");
        retVal.in("Hammer.hammer", size);
        return retVal;
    }

    /**
     * @return a batch hammer query for presence only
     *
     * @param size	number of kmers to be passed in
     *
     * @throws SQLException
     */
    protected DbQuery buildBatchSetQuery(int size) throws SQLException {
        DbQuery retVal = new DbQuery(this.db, "Hammer");
        retVal.select("Hammer", "hammer");
        retVal.in("Hammer.hammer", size);
        return retVal;
    }

    /**
     * Run a database query to find all the genomes for the specified kmers and update the score map.
     *
     * @param query		query to fill with parameters and run
     * @param kmers		parameter values to store in the query, indicating the kmers to run, mapped
     * 					to the number of hits for each
     * @param gCounts	score map for recording results
     *
     * @throws SQLException
     */
    private void runCountQuery(DbQuery query, CountMap<String> kmers, WeightMap gCounts) throws SQLException {
        this.setupBatchQuery(query, kmers.keys());
        int count = 0;
        for (DbRecord result : query) {
            HammerDb.Source source = new HammerDb.Source(result.getString("Hammer.fid"),
                    result.getString("Hammer.role"), result.getDouble("Hammer.strength"));
            String hammer = result.getString("Hammer.hammer");
            this.countHit(gCounts, source, kmers.getCount(hammer));
            count++;
        }
        log.debug("{} kmers queried and {} results found.", kmers.size(), count);
    }

    /**
     * Set up a batch query to ask for the specified kmers.
     *
     * @param query		query to set up
     * @param kmers		set of kmers to query
     *
     * @throws SQLException
     */
    protected void setupBatchQuery(DbQuery query, Set<String> kmers) throws SQLException {
        int idx = 1;
        for (String kmer : kmers) {
            query.setParm(idx, kmer);
            idx++;
        }
    }

    @Override
    protected void findHitsInternal(Collection<Hit> collection, Collection<? extends ISequence> seqs, int kSize, boolean dir) {
        // Again, we batch the queries for performance.  A batch consists of a map from each kmer to a list of hits
        // with dummy fids and roles. The output collection will retain the hits that are found.
        try {
            // Set up a map to hold the batches.
            Map<String, List<HammerDb.Hit>> batchMap = new HashMap<String, List<HammerDb.Hit>>(this.hashSize);
            // Create a query.
            DbQuery query = null;
            try {
                query = this.buildBatchQuery(this.batchSize);
                // We loop through the sequences, then loop through the kmers in each sequence.  We can't use any
                // of the kmer iterators, because we need to know where each hammer was found.
                for (ISequence seq : seqs) {
                    String seqId = seq.getLabel();
                    String dna = seq.getSequence();
                    final int len = dna.length();
                    final int n = len - kSize;
                    for (int i = 0; i < n; i++) {
                        if (batchMap.size() == this.batchSize) {
                            // Now we run the query.
                            this.runHitQuery(query, batchMap, collection);
                            // Clear the kmer set for the next query.
                            batchMap.clear();
                        }
                        // Save this kmer in the map.
                        var hit = new HammerDb.Hit(seqId, len, i, dir, "", "", kSize, 1.0);
                        List<HammerDb.Hit> hitList = batchMap.computeIfAbsent(dna.substring(i, i + kSize),
                                x -> new ArrayList<HammerDb.Hit>(5));
                        hitList.add(hit);
                    }
                }
            } finally {
                if (query != null)
                    query.close();
            }
            // Do we have a residual?
            if (! batchMap.isEmpty()) {
                // Yes.  Create a mini-query for it.
                query = null;
                try {
                    query = this.buildBatchQuery(batchMap.size());
                    this.runHitQuery(query, batchMap, collection);
                } finally {
                    if (query != null)
                        query.close();
                }
            }

        } catch (SQLException e) {
            throw new RuntimeException("Find hammer hits error:  " + e.getMessage());
        }
    }

    /**
     * Run a query to test a set of kmers for presence in the hammer database.  If any are found, the feature
     * IDs will be filled into the hit objects and they will be added to the output.
     *
     * @param query			query object to use
     * @param batchMap		map of kmers to lists of hit objects
     * @param collection	output collection in which to store the hits
     *
     * @throws SQLException
     */
    private void runHitQuery(DbQuery query, Map<String, List<HammerDb.Hit>> batchMap, Collection<HammerDb.Hit> collection) throws SQLException {
        this.setupBatchQuery(query, batchMap.keySet());
        int count = 0;
        for (DbRecord result : query) {
            String fid = result.getString("Hammer.fid");
            double strength = result.getDouble("Hammer.strength");
            String role = result.getString("Hammer.role");
            // Now we get the hits, update the fids, and add the hits to the output.
            List<HammerDb.Hit> hitList = batchMap.get(result.getString("Hammer.hammer"));
            hitList.stream().forEach(x -> x.setHammerSource(fid, role, strength));
            collection.addAll(hitList);
            count++;
        }
        log.info("{} kmers queried, {} were hammers.", batchMap.size(), count);
    }

    @Override
    public File getLoadFile() {
        // The load file is not available for this type of hammer database.
        return null;
    }

    @Override
    protected void findHammersInternal(HashSet<String> hammerSet, String seq, int kSize) {
        // For performance reasons, we batch the kmer queries.  We will collect a batch of
        // kmers, ask for the feature IDs, and count the results.  This gets repeated
        // until we reach the end.
        try {
            // Get an iterator for the kmers.
            SequenceKmerIterable kIter = new SequenceKmerIterable(seq, kSize);
            // This will hold the current kmer batch.
            var kmers = new HashSet<String>(this.hashSize);
            // Create a query.
            DbQuery query = null;
            try {
                 query = this.buildBatchSetQuery(this.batchSize);
                // Loop through the kmers, building a batch.
                for (String kmer : kIter) {
                    kmers.add(kmer);
                    if (kmers.size() == this.batchSize) {
                        // Now we run the query.
                        this.runHammerQuery(query, kmers, hammerSet);
                        // Clear the kmer set for the next query.
                        kmers.clear();
                    }
                }
            } finally {
                if (query != null)
                    query.close();
            }
            // Do we have a residual?
            if (kmers.size() > 0) {
                // Yes.  Create a mini-query for it.
                query = null;
                try {
                    query = this.buildBatchSetQuery(kmers.size());
                    this.runHammerQuery(query, kmers, hammerSet);
                } finally {
                    if (query != null)
                        query.close();
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException("Find-closest query error: " + e.toString());
        }

    }

    /**
     * Run a query to determine which kmers in a set are hammers.
     *
     * @param query		query object to use
     * @param kmers		set of kmers to check
     * @param hammerSet	output set into which hammers found are placed
     *
     * @throws SQLException
     */
    private void runHammerQuery(DbQuery query, HashSet<String> kmers, HashSet<String> hammerSet) throws SQLException {
        this.setupBatchQuery(query, kmers);
        int count = 0;
        for (DbRecord result : query) {
            String hammer = result.getString("Hammer.hammer");
            kmers.add(hammer);
        }
        log.info("{} kmers queried, {} were hammers.", kmers.size(), count);

    }

    @Override
    public Map<String, Source> findGenomeHammers(String genomeId) {
        // Set up the output hash.
        Map<String, Source> retVal = new HashMap<String, Source>();
        // Create a query to look for all hammers with a matching feature ID.
        try (DbQuery query = new DbQuery(this.db, "Hammer")) {
            query.select("Hammer", "fid", "role", "hammer", "strength");
            query.rel("Hammer.fid", Relop.LIKE);
            query.setParm(1, "fig|" + genomeId + ".%");
            // Loop through the results, adding them to the return hash.
            for (DbRecord result : query) {
                String hammer = result.getString("Hammer.hammer");
                Source source = new Source(result.getString("Hammer.fid"), result.getString("Hammer.role"),
                        result.getDouble("Hammer.strength"));
                retVal.put(hammer, source);
            }
        } catch (SQLException e) {
            // Rethrow as a run-time exception.
            throw new RuntimeException(e);
        }
        return retVal;
    }

}
