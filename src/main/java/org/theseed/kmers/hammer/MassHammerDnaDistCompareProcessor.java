/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.io.TabbedLineReader;
import org.theseed.kmers.filter.CompareFilter;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.proteins.hammer.HashHammerDb;
import org.theseed.reports.HammerDnaDistReporter;
import org.theseed.utils.BaseMultiReportProcessor;

/**
 * This is a utility command that runs the HammerDnaDistCompareProcessor on multiple input sets. Each input
 * set is run with one or more user-specified report formats and the results placed in the specified
 * output directory. The hammers will be loaded once, removing the biggest time sink in the individual
 * reporting tasks.
 *
 * The standard input should contain a list of input file names to use. The input should be tab-delimited,
 * with headers, the input file names being in the first line. The output file names will
 * be based on the input file base name with the extension removed and the report type integrated.
 *
 * The positional parameters should be the name of the hammer database file and the name of the finder
 * directory used to create the hammers. The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -D	output directory name (default "compare" in the current directory)
 * -i	input file (if not STDIN)
 *
 * --clear		erase the output directory before processing
 * --filter		type of comparison filtering to use (default EXTREME)
 * --format		type of report to produce (multiple, default all of them)
 *
 *
 * @author Bruce Parrello
 *
 */
public class MassHammerDnaDistCompareProcessor extends BaseMultiReportProcessor {

	// FIELDS
	/** logging facility */
	private static final Logger log = LoggerFactory.getLogger(MassHammerDnaDistCompareProcessor.class);
	/** input file set */
	private Set<File> inFiles;
	/** command processor to produce the reports */
	private HammerDnaDistCompareProcessor subProcessor;
	/** master hammer database */
	private HammerDb hammers;

	// COMMAND-LINE OPTIONS

	/** type of filtering to use */
	@Option(name = "--filter", usage = "type of comparison filter to use")
	private CompareFilter.Type filterType;

	/** list of report formats */
	@Option(name = "--format", usage = "type of desired output report (multiple)")
	private List<HammerDnaDistReporter.Type> reportTypes;

	/** input file (if not STDIN) */
	@Option(name = "--input", aliases = { "-i" }, metaVar = "inFiles.tbl",
			usage = "file containing input file names (if not STDIN)")
	private File inFile;

	/** name of the master hammer file */
	@Argument(index = 0, metaVar = "hammers.tbl", usage = "name of the hammer database file", required = true)
	private File hammerFile;

	/** name of the protein finder directory used to create the hammers */
	@Argument(index = 1, metaVar = "finderDir", usage = "name of the protein finder directory", required = true)
	private File finderDir;

	@Override
	protected File setDefaultOutputDir(File curDir) {
		return new File(curDir, "compare");
	}

	@Override
	protected void setMultiReportDefaults() {
		// Default to no reports.
		this.reportTypes = new ArrayList<HammerDnaDistReporter.Type>();
		// Use extreme-style filtering.
		this.filterType = CompareFilter.Type.EXTREME;
		// Default to STDIN.
		this.inFile = null;
	}

	@Override
	protected void validateMultiReportParms() throws IOException, ParseFailureException {
		// If no reports were specified, specify all of them.
		if (this.reportTypes.isEmpty())
			Arrays.stream(HammerDnaDistReporter.Type.values()).forEach(x -> this.reportTypes.add(x));
		// First, we must verify the input file list.
		this.inFiles = new TreeSet<File>();
		try (TabbedLineReader inStream = this.openInput()) {
			for (var line : inStream) {
				File inFile1 = new File(line.get(0));
				if (! inFile1.canRead())
					throw new FileNotFoundException("Proposed input file " + inFile1.getAbsolutePath()
							+ " is not found or unreadable.");
				this.inFiles.add(inFile1);
			}
		}
		log.info("{} input files found.", this.inFiles.size());
		// Now verify the finder directory.
		if (! this.finderDir.isDirectory())
			throw new FileNotFoundException("Finder directory " + this.finderDir + " is not found or invalid.");
		// Initialize the common parameters.
		this.subProcessor = new HammerDnaDistCompareProcessor();
		this.subProcessor.setupAllRuns(this.finderDir, this.filterType);
		// Now load the hammer database.
		if (! this.hammerFile.canRead())
			throw new FileNotFoundException("Hammer file " + this.hammerFile + " is not found or unreadable.");
		this.hammers = new HashHammerDb(this.hammerFile);
		// Finally, inform the user of the reports selected.
		log.info("Reports selected: {}.", StringUtils.join(this.reportTypes, ", "));
	}

	/**
	 * Open the input stream.
	 *
	 * @return the open input stream as a tabbed line reader
	 *
	 * @throws IOException
	 */
	private TabbedLineReader openInput() throws IOException {
		TabbedLineReader retVal;
		if (this.inFile == null) {
			log.info("Reading input files from the standard input.");
			retVal = new TabbedLineReader(System.in);
		} else if (! this.inFile.canRead())
			throw new FileNotFoundException("Input list file " + this.inFile + " is not found or unreadable.");
		else {
			log.info("Reading input file list from {}.", this.inFile);
			retVal = new TabbedLineReader(this.inFile);
		}
		return retVal;
	}

	@Override
	protected void runMultiReports() throws Exception {
		// This will count the number of input files processed.
		int fileCount = 0;
		final int fileTotal = this.inFiles.size();
		// This will count the number of reports written.
		int outCount = 0;
		// Loop through the input files.
		for (File inFile1 : this.inFiles) {
			// Load the input file.
			fileCount++;
			log.info("Processing input file {} of {}: {}.", fileCount, fileTotal, inFile1);
			Set<String> genomes = TabbedLineReader.readSet(inFile1, "1");
			// Loop through the report types.
			for (var reportType1 : this.reportTypes) {
				// Set up the report.
				this.subProcessor.setupOneRun(genomes, reportType1);
				// Now we calculate the output file name.
				String inBase = StringUtils.substringBeforeLast(inFile1.getName(), ".");
				String outName = inBase + "." + StringUtils.lowerCase(reportType1.toString()) + ".tbl";
				File outFile = this.getOutFile(outName);
				log.info("Producing {} report for {} on {}.", reportType1, inFile1, outFile);
				// Open the output file and generate the report.
				try (PrintWriter writer = new PrintWriter(outFile)) {
					this.subProcessor.runHammers(hammers, writer);
					outCount++;
				}
			}
		}
		log.info("{} reports produced for {} input files.", outCount, fileCount);
	}

}
