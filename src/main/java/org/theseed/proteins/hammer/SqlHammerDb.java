/**
 *
 */
package org.theseed.proteins.hammer;

import java.io.File;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.theseed.counters.CountMap;
import org.theseed.java.erdb.DbConnection;
import org.theseed.java.erdb.DbLoader;
import org.theseed.java.erdb.DbQuery;
import org.theseed.java.erdb.DbRecord;
import org.theseed.java.erdb.SqlBuffer;
import org.theseed.sequence.Sequence;

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
                throw new RuntimeException("No hammers found in hammer database.");
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
        public void updateHammerMap(String genome, String hammer) {
            try {
                this.loader.set("hammer", hammer);
                this.loader.set("genome_id", genome);
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
    protected CountMap<String> findClosestInternal(Collection<Sequence> seqs, int kSize) {
        CountMap<String> retVal = new CountMap<String>();
        // For performance reason, we batch the kmer queries.  We will collect a batch of
        // kmers, ask for the genome IDs, and count the results.  This gets repeated
        // until we reach the end.
        try {
            // Get an iterator for the kmers.
            Iterable<String> kIter = KmerSeries.init(seqs, kSize);
            log.debug("Iterating through {} sequences.", seqs.size());
            // This will hold the current kmer batch.
            Set<String> kmers = new HashSet<String>(this.hashSize);
            // Create a query.
            try (DbQuery query = new DbQuery(this.db, "Hammer")) {
                query.select("Hammer", "genome_id");
                query.in("Hammer.hammer", this.batchSize);
                // Loop through the kmers, building a batch.
                for (String kmer : kIter) {
                    kmers.add(kmer);
                    if (kmers.size() == this.batchSize) {
                        // Now we run the query.
                        this.runQuery(query, kmers, retVal);
                        // Clear the kmer set for the next query.
                        kmers.clear();
                    }
                }
            }
            // Do we have a residual?
            if (! kmers.isEmpty()) {
                // Yes.  Create a mini-query for it.
                try (DbQuery query = new DbQuery(this.db, "Hammer")) {
                    query.select("Hammer", "genome_id");
                    query.in("Hammer.hammer", kmers.size());
                    this.runQuery(query, kmers, retVal);
                }
            }
        } catch (SQLException e) {
            throw new RuntimeException("Find-closest query error: " + e.toString());
        }
        return retVal;
    }

    /**
     * Run a database query to find all the genomes for the specified kmers and update the count map.
     *
     * @param query		query to fill with parameters and run
     * @param kmers		parameter values to store in the query, indicating the kmers to run
     * @param gCounts	count map for recording results
     *
     * @throws SQLException
     */
    private void runQuery(DbQuery query, Set<String> kmers, CountMap<String> gCounts) throws SQLException {
        int idx = 1;
        for (String kmer : kmers) {
            query.setParm(idx, kmer);
            idx++;
        }
        int count = 0;
        for (DbRecord result : query) {
            gCounts.count(result.getString("Hammer.genome_id"));
            count++;
        }
        log.debug("{} kmers queried and {} results found.", idx-1, count);
    }

}
