/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.WeightMap;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command produces a bin report based on Kraken output files in a sample group.  The Kraken files should be
 * in the sample subdirectories with the name "kraken.tbl".  The input directory should be a master sample directory;
 * that is, a directory whose subdirectories each contain sample data.  If the desired data file is not found, the
 * subdirectory is simply skipped.
 *
 * The positional parameters are the name of the input master directory, the name of the repgen taxonomy weight
 * file for the repgen set used to build the hammers, and the repgen stats file for the same set.  The taxonomy
 * weights allow us to convert the kraken taxon IDs to repgen IDs, and the stats file gives us the name of each
 * representative genome.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for bin report (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class KrakenBinReportProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(KrakenBinReportProcessor.class);
    /** list of sample directories */
    private File[] sampleDirs;
    /** map of repgen IDs to repgen names */
    private Map<String, String> nameMap;
    /** map from genus ID -> repgen weight map */
    private Map<Integer, WeightMap> genusMaps;
    /** map from species ID -> repgen weight map */
    private Map<Integer, WeightMap> speciesMaps;
    /** file filter for sample directories */
    private FileFilter SAMPLE_DIR_FILTER = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            boolean retVal = pathname.isDirectory();
            if (retVal) {
                File reportFile = new File(pathname, "kraken.tbl");
                retVal = reportFile.isFile();
            }
            return retVal;
        }
    };

    // COMMAND-LINE OPTIONS

    /** name of the master sample directory */
    @Argument(index = 0, metaVar = "inDir", usage = "master sample directory", required = true)
    private File inDir;

    /** name of the repgen taxonomic weight file */
    @Argument(index = 1, metaVar = "repXX.tax.tbl", usage = "name of the taxonomic weight file for the appropriate repgen set",
            required = true)
    private File repTaxFile;

    /** name of the repgen stats file */
    @Argument(index = 2, metaVar = "repXX.stats.tbl", usage = "name of the statistics file for the appropriate repgen set",
            required = true)
    private File repStatsFile;


    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Get the sample directories from the input master.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input master directory " + this.inDir + " is not found or invalid.");
        this.sampleDirs = this.inDir.listFiles(SAMPLE_DIR_FILTER);
        if (this.sampleDirs.length == 0)
            throw new IOException("Input master directory " + this.inDir + " does not contain any Kraken reports.");
        log.info("{} samples found with Kraken reports in {}.", this.sampleDirs.length, this.inDir);
        // Verify the two data files.
        if (! this.repTaxFile.canRead())
            throw new FileNotFoundException("Taxonomic weight file " + this.repTaxFile + " is not found or unreadable.");
        if (! this.repStatsFile.canRead())
            throw new FileNotFoundException("Statistics file " + this.repStatsFile + " is not found or unreadable.");
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // First, read in the genome names.
        this.nameMap = TabbedLineReader.readMap(this.repStatsFile, "rep_id", "rep_name");
        log.info("{} genome names read from {}.", this.nameMap.size(), this.repStatsFile);
        // Second, read in the taxonomic weights.
        this.readTaxFile();
        // Compute the optimal hash size.
        int hashSize = this.nameMap.size() * 4 / 3 + 1;
        // Write the report header.
        writer.println("sample_id\trepgen_id\trepgen_name\tcount");
        // Now we loop through the samples, using the bin report to compute weights.
        for (File sampleDir : this.sampleDirs) {
            String sampleId = sampleDir.getName();
            log.info("Processing sample {}.", sampleId);
            WeightMap weights = new WeightMap(hashSize);
            // Loop through the Kraken output file.  This is a tab-delimited file with no headers.
            final File krakenFile = new File(sampleDir, "kraken.tbl");
            try (TabbedLineReader kStream = new TabbedLineReader(krakenFile, 5)) {
                int lineCount = 0;
                int genusCount = 0;
                int genusNotFound = 0;
                int speciesCount = 0;
                int speciesNotFound = 0;
                for (var line : kStream) {
                    int taxId = line.getInt(4);
                    lineCount++;
                    // Get the rank for this taxonomic grouping.
                    String rankCode = line.get(3);
                    switch (rankCode) {
                    case "G" :
                        // For a genus, we count only the taxon-specific fragments.
                        WeightMap gWeights = this.genusMaps.get(taxId);
                        if (gWeights == null)
                            genusNotFound++;
                        else {
                            weights.accumulate(gWeights, line.getInt(2));
                            genusCount++;
                        }
                        break;
                    case "S" :
                        // For a species, we cound the whole clade, since we don't go any lower.
                        WeightMap sWeights = this.speciesMaps.get(taxId);
                        if (sWeights == null)
                            speciesNotFound++;
                        else {
                            weights.accumulate(sWeights, line.getInt(1));
                            speciesCount++;
                        }
                        break;
                    }
                }
                log.info("{} lines read, {} species not found, {} genera not found, {} species counted, {} genera counted, {} representatives.",
                        lineCount, speciesNotFound, genusNotFound, speciesCount, genusCount, weights.size());
            }
            log.info("Writing scores for {}.", sampleId);
            var scores = weights.sortedCounts();
            for (var score : scores) {
                String repId = score.getKey();
                double count = score.getCount();
                String repName = this.nameMap.getOrDefault(repId, "<unknown>");
                writer.println(sampleId + "\t" + repId + "\t" + repName + "\t" + count);
            }
        }

    }

    /**
     * Read in the taxonomic weights file and build a weight map for each genus and species.
     *
     * @throws IOException
     */
    private void readTaxFile() throws IOException {
        log.info("Loading taxonomic weights from {}.", this.repTaxFile);
        // Create the blank maps.
        this.genusMaps = new HashMap<Integer, WeightMap>(8000);
        this.speciesMaps = new HashMap<Integer, WeightMap>(50000);
        try (TabbedLineReader  taxStream = new TabbedLineReader(this.repTaxFile)) {
            int taxColIdx = taxStream.findField("tax_id");
            int rankColIdx = taxStream.findField("rank");
            int repColIdx = taxStream.findField("repgen_id");
            int weightColIdx = taxStream.findField("weight");
            // Loop through the taxonomy file, building weight maps.
            for (var line : taxStream) {
                // Get the main fields.
                int taxId = line.getInt(taxColIdx);
                String repgenId = line.get(repColIdx);
                double weight = line.getDouble(weightColIdx);
                // Compute the target map based on rank.
                Map<Integer, WeightMap> targetMap;
                if (line.get(rankColIdx).equals("genus"))
                    targetMap = this.genusMaps;
                else
                    targetMap = this.speciesMaps;
                // Get the weight map for the specified taxon ID.
                WeightMap weights = targetMap.computeIfAbsent(taxId, x -> new WeightMap());
                // Add the weight.
                weights.count(repgenId, weight);
            }
        }
        log.info("{} species and {} genera tabulated.", this.speciesMaps.size(), this.genusMaps.size());
    }

}
