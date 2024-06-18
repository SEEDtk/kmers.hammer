/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseReportProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;

/**
 * This command reads multiple files output from the ReadTestProcessor and outputs various statistics.  The idea here is that
 * we expect a file of good hits and a file of bad hits, and we want to compare basic statistics for both.  The output will
 * list the number of hits for each role, the weighted number of hits for each role, and the minimum, maximum, mean, and
 * standard deviation for the weights and qualities.
 *
 * The positional parameters are the names of the input files.  All the files must have the same headers in the same order,
 * and these need to include "role", "weight", and "quality".
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file for report (if not STDOUT)
 *
 * @author Bruce Parrello
 *
 */
public class HitStatsProcessor extends BaseReportProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(HitStatsProcessor.class);
    /** map of stream names to input streams */
    private SortedMap<String, TabbedLineReader> inStreams;
    /** map of stream names to role-ids to hit-weight statistics */
    private SortedMap<String, Map<String, SummaryStatistics>> roleMaps;
    /** map of stream names to total weight statistics */
    private SortedMap<String, SummaryStatistics> weightMap;
    /** map of stream names to quality statistics */
    private SortedMap<String, SummaryStatistics> qualMap;
    /** column index of the weight column */
    private int weightColIdx;
    /** column index of the quality column */
    private int qualColIdx;
    /** column index of the role-id column */
    private int roleColIdx;
    /** name of the first input stream */
    private String firstName;
    /** null summary-statistics object */
    private static final SummaryStatistics NULL_STATS = new SummaryStatistics();

    // COMMAND-LINE OPTIONS

    /** list of input files */
    @Argument(index = 0, metaVar = "inFile1 inFile2 ...", usage = "input files containing hammer hit data", required = true)
    private List<File> inFiles;

    @Override
    protected void setReporterDefaults() {
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Loop through the input files and verify that they are readable.  We also verify that each one has a unique
        // base name.
        var fileMap = new TreeMap<String, File>();
        for (File inFile : this.inFiles) {
            if (! inFile.canRead())
                throw new FileNotFoundException("Input file " + inFile + " is not found or invalid.");
            // Compute the base name for this file.
            String baseName = inFile.getName();
            if (fileMap.containsKey(baseName))
                throw new ParseFailureException("Input file " + inFile + " has the same base name as another input file.");
            fileMap.put(baseName, inFile);
        }
        // Now that we know the input files are all valid, we can open them and put them in the maps.
        // Create the four maps.
        this.inStreams = new TreeMap<String, TabbedLineReader>();
        this.roleMaps = new TreeMap<String, Map<String, SummaryStatistics>>();
        this.qualMap = new TreeMap<String, SummaryStatistics>();
        this.weightMap = new TreeMap<String, SummaryStatistics>();
        // Protect us from errors.
        try {
            this.firstName = null;
            // Now loop through the files.
            for (var fileEntry : fileMap.entrySet()) {
                // Get the base name and open the associated file.
                String baseName = fileEntry.getKey();
                TabbedLineReader fileStream = new TabbedLineReader(fileEntry.getValue());
                this.inStreams.put(baseName, fileStream);
                // Is this the first stream?
                if (firstName == null) {
                    this.firstName = baseName;
                    // Here it is the first stream, so save the column indices.
                    this.qualColIdx = fileStream.findField("quality");
                    this.weightColIdx = fileStream.findField("weight");
                    this.roleColIdx = fileStream.findField("role");
                } else {
                    // Here it is a subsequent stream, so verify the column indices.
                    if (fileStream.findColumn("quality") != this.qualColIdx)
                        throw new IOException("Quality field must be in column " + this.qualColIdx + " in " + baseName + " file.");
                    if (fileStream.findColumn("weight") != this.weightColIdx)
                        throw new IOException("Weight field must be in column " + this.weightColIdx + " in " + baseName + " file.");
                    if (fileStream.findColumn("role") != this.roleColIdx)
                        throw new IOException("Role-id field must be in column " + this.roleColIdx + " in " + baseName + " file.");
                }
                // Create the statistical receptacles for this stream.
                this.roleMaps.put(baseName, new HashMap<String, SummaryStatistics>(30));
                this.qualMap.put(baseName, new SummaryStatistics());
                this.weightMap.put(baseName, new SummaryStatistics());
            }
            log.info("{} input files prepared for processing.", this.inStreams.size());
        } catch (IOException e) {
            this.close();
            throw new IOException(e);
        }
    }

    /**
     *
     */
    protected void close() {
        // Close all the open streams.
        for (TabbedLineReader inStream : this.inStreams.values())
            inStream.close();
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        try {
            // Loop through the input files.  Note that because we have pre-filled all the maps
            // at stream level, we can do this in parallel.
            log.info("Processing input streams.");
            this.inStreams.entrySet().stream().forEach(x -> this.analyzeFile(x));
            log.info("Producing output report.");
            // We will buffer our output lines in here.
            StringBuilder outBuffer = new StringBuilder(80);
            // Build the header line.
            outBuffer.append("value");
            for (String baseName : inStreams.keySet())
                outBuffer.append('\t').append(baseName);
            writer.println(outBuffer.toString());
            // Write out the quality and weight data.
            this.writeNumData(writer, "quality", this.qualMap);
            this.writeNumData(writer, "weight", this.weightMap);
            // Now, invert the rolemaps so for each role we have a map from stream name to stats.  First, we
            // need a comprehensive list of roles IDs.
            Set<String> roles = new TreeSet<String>();
            for (var roleMap : this.roleMaps.values())
                roles.addAll(roleMap.keySet());
             // Now we loop through the roles, constructing maps to send to the output writer.
            for (String role : roles) {
                // Construct the map for this role.
                SortedMap<String, SummaryStatistics> roleMap = new TreeMap<String, SummaryStatistics>();
                for (var streamEntry : this.roleMaps.entrySet()) {
                    String baseName = streamEntry.getKey();
                    roleMap.put(baseName, streamEntry.getValue().getOrDefault(role, NULL_STATS));
                }
                this.writeNumData(writer, role, roleMap);
            }
            // Insure all the output has been written.
            writer.flush();
            log.info("Report complete.");
        } finally {
            // Close all the open streams.
            this.close();
        }
    }

    /**
     * Process all the records in the specified input file and accumulate its statistics.
     *
     * @param fileEntry		stream map entry containing the stream name and the stream itself
     */
    private void analyzeFile(Entry<String, TabbedLineReader> fileEntry) {
        // Parse the incoming stream map entry.
        String baseName = fileEntry.getKey();
        TabbedLineReader inStream = fileEntry.getValue();
        // Get the statistic objects and the role map.
        SummaryStatistics qualStats = this.qualMap.get(baseName);
        SummaryStatistics weightStats = this.weightMap.get(baseName);
        Map<String, SummaryStatistics> roleMap = this.roleMaps.get(baseName);
        // Let the user know we're reading this file.
        log.info("Processing input file {}.", baseName);
        long lastMsg = System.currentTimeMillis();
        // Loop through the input lines.
        for (var line : inStream) {
            String roleId = line.get(this.roleColIdx);
            double weight = line.getDouble(this.weightColIdx);
            double quality = line.getDouble(this.qualColIdx);
            // Record the role hit.
            SummaryStatistics roleStats = roleMap.computeIfAbsent(roleId, x -> new SummaryStatistics());
            roleStats.addValue(weight);
            // Record the file-wide hit.
            qualStats.addValue(quality);
            weightStats.addValue(weight);
            long now = System.currentTimeMillis();
            if (now - lastMsg >= 10000) {
                log.info("{} records read in {}.", weightStats.getN(), baseName);
                lastMsg = now;
            }
        }
        log.info("{} total records read in {}.", weightStats.getN(), baseName);
    }

    /**
     * Write out a statistical analysis for a specified quantity.
     *
     * @param writer		output print writer
     * @param label			label of the quantity
     * @param dataMap		map of stream names to summary-statistics objects
     */
    private void writeNumData(PrintWriter writer, String label, SortedMap<String, SummaryStatistics> dataMap) {
        // Write the count.
        writer.println(dataMap.values().stream().map(x -> Long.toString(x.getN())).collect(Collectors.joining("\t", label + " count", "")));
        // Write the minimum.
        writer.println(dataMap.values().stream().map(x -> Double.toString(x.getMin())).collect(Collectors.joining("\t", label + " min", "")));
        // Write the mean.
        writer.println(dataMap.values().stream().map(x -> Double.toString(x.getMean())).collect(Collectors.joining("\t", label + " mean", "")));
        // Write the maximum.
        writer.println(dataMap.values().stream().map(x -> Double.toString(x.getMax())).collect(Collectors.joining("\t", label + " max", "")));
        // Write the standard deviation.
        writer.println(dataMap.values().stream().map(x -> Double.toString(x.getStandardDeviation())).collect(Collectors.joining("\t", label + " sdev", "")));
    }

}
