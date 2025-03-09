/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.kmers.filter.CompareFilter;
import org.theseed.proteins.hammer.HammerDb;
import org.theseed.reports.HammerDnaDistReporter;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.seeds.ProteinFinder;
import org.theseed.utils.BaseHammerUsageProcessor;

/**
 * This command takes in a set of genomes and compares the DNA kmer distance to the hammer distance
 * for each of the SOURs.
 *
 * The comparison is performed using the finder that generated the hammer set. A list of genome IDs
 * must be presented on the standard input. We then read in the finder file for each SOUR and keep
 * the DNA sequences for each genome in the input list.  These sequences are then compared for hammer
 * distance and DNA kmer distance (the kmer size being the same as the hammer size). The default
 * output report shows various correlation coefficients between the two distance types for each SOUR.
 *
 * The positional parameter should be the name of the finder directory used to create the hammers.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	output file (if not STDOUT)
 * -b	batch size for queries
 * -i	input file containing the genome IDs (if not STDIN)
 * -c	index (1-based) or name of the input column containing the genome IDs (default "1")
 *
 * --hType		type of hammer database (default MEMORY)
 * --method		voting method to use (default COUNT)
 * --file		file containing hammer database (either SQLite database or hammer flat file)
 * --url		URL of database (host and name, MySQL only)
 * --parms		database connection parameter string (MySQL only)
 * --type		database engine type
 * --format		format of the output report (default CORRELATION)
 * --filter		filtering method to use for comparisons (default EXTREME)
 *
 * @author Bruce Parrello
 *
 */
public class HammerDnaDistCompareProcessor extends BaseHammerUsageProcessor
		implements HammerDnaDistReporter.IParms, CompareFilter.IParms {

	// FIELDS
	/** logging facility */
	protected static Logger log = LoggerFactory.getLogger(HammerDnaDistCompareProcessor.class);
	/** output report writer */
	private HammerDnaDistReporter reporter;
	/** protein finder */
	private ProteinFinder finder;
	/** set of genome IDs */
	private Set<String> genomeSet;
	/** trivial-compare filter */
	private CompareFilter filter;
	/** map of SOURs to FASTA files */
	private Map<String, File> sourMap;

	// COMMAND-LINE OPTIONS

	/** input file containing genome IDs */
	@Option(name = "--input", aliases = { "-i" }, metaVar = "genomeIDs.tbl",
			usage = "name of the input file containing the genome IDs (if not STDIN)")
	private File inFile;

	/** index (1-based) or name of input column containing genome IDs */
	@Option(name = "--col", aliases = { "-c" }, metaVar = "id_col",
			usage = "index (1-based) or name of the column containing the genome IDs")
	private String idCol;

	/** format for the output report */
	@Option(name = "--format", usage = "format of the output report")
	private HammerDnaDistReporter.Type reportType;

	/** type of comparison filtering to use */
	@Option(name = "--filter", usage = "method for filtering comparisons")
	private CompareFilter.Type filterType;

	/** name of the finder directory used to create the hammers */
	@Argument(index = 0, metaVar = "finderDir", usage = "name of the finder directory used to create the hammers",
			required = true)
	private File finderDir;

	@Override
	protected void setHammerDefaults() {
		this.inFile = null;
		this.idCol = "1";
		this.reportType = HammerDnaDistReporter.Type.CORRELATION;
		this.filterType = CompareFilter.Type.EXTREME;
	}

	@Override
	protected void validateHammerParms() throws IOException, ParseFailureException {
		// Load the protein finder.
		if (! this.finderDir.isDirectory())
			throw new FileNotFoundException("Finder directory " + this.finderDir + " is not found or invalid.");
		this.finder = new ProteinFinder(this.finderDir);
		// Create the comparison filter. We must do this before creating the report writer.
		this.filter = this.filterType.create(this);
		// Now read in the genome IDs.
		this.genomeSet = new HashSet<String>();
		try (TabbedLineReader gReader = this.openInput()) {
			// Get the index of the genome ID column.
			int idColIdx = gReader.findField(this.idCol);
			for (var line : gReader) {
				String genomeId = line.get(idColIdx);
				this.genomeSet.add(genomeId);
			}
		}
		log.info("{} genome IDs read from input.", this.genomeSet.size());
		// Create the output report writer.
		this.reporter = this.reportType.create(this);
	}

	/**
	 * Set up the parameters for all externally-commanded runs.
	 *
	 * @param finderDir1	the protein-finder directory to use
	 * @param filterType1	the type of comparison filtering to use
	 *
	 * @throws IOException
	 */
	public void setupAllRuns(File finderDir1, CompareFilter.Type filterType1) throws IOException {
		this.finder = new ProteinFinder(finderDir1);
		log.info("Using finder in directory {}.", finderDir1);
		this.filter = filterType1.create(this);
		log.info("Filtering out {} comparisons.", filterType1);
		// Get the map of SOURs to FASTA files.
		this.sourMap = this.finder.getFastas();
		log.info("{} SOURs in protein finder.", sourMap.size());
	}

	/**
	 * Set up the parameters for an individual run.
	 *
	 * @param gSet			input genome set to use
	 * @param reportType1	type of report to write
	 */
	public void setupOneRun(Set<String> gSet, HammerDnaDistReporter.Type reportType1) {
		this.reporter = reportType1.create(this);
		this.genomeSet = gSet;
		log.info("Producing {} report for {}-genome input set.", reportType1, gSet.size());
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
			log.info("Reading genome IDs from the standard input.");
			retVal = new TabbedLineReader(System.in);
		} else if (! this.inFile.canRead())
			throw new FileNotFoundException("Genome ID file " + this.inFile + " is not found or unreadable.");
		else {
			log.info("Reading genome IDs from {}.", this.inFile);
			retVal = new TabbedLineReader(this.inFile);
		}
		return retVal;
	}

	@Override
	protected void runHammers(HammerDb hammerDb, PrintWriter writer) throws Exception {
		// Start the report.
		this.reporter.openReport(writer);
		// Save the kmer size for performance.
		DnaKmers.setKmerSize(hammerDb.getKmerSize());
		// Loop through the SOURs.
		for (var sourEntry : this.sourMap.entrySet()) {
			// Get the data for this SOUR.
			String sourName = sourEntry.getKey();
			File sourFile = sourEntry.getValue();
			// Set up the report.
			this.reporter.openSourReport(sourName);
			log.info("Scanning {} for {} proteins.", sourFile, sourName);
			// We'll accumulate the desired sequence kmers in here.
			NavigableMap<String, DnaKmers> dnaMap = new TreeMap<String, DnaKmers>();
			// We'll accumulate the hammer sets in here.
			Map<String, Set<String>> hammerMap = new HashMap<String, Set<String>>(this.genomeSet.size() * 5 / 3 + 1);
			// Now extract the sequences for the desired genomes and put their kmers in the map.
			try (FastaInputStream seqStream = new FastaInputStream(sourFile)) {
				for (Sequence seq : seqStream) {
					// Compute the genome ID.
					String genomeID = Feature.genomeOf(seq.getLabel());
					if (this.genomeSet.contains(genomeID)) {
						// In rare cases we might have two sequences for the same genome. When that happens,
						// we keep the longest.
						DnaKmers oldKmers = dnaMap.get(genomeID);
						if (oldKmers == null || oldKmers.getDna().length() < seq.getSequence().length()) {
							DnaKmers newKmers = new DnaKmers(seq.getSequence());
							dnaMap.put(genomeID, newKmers);
							Set<String> hammerSet = hammerDb.findHammers(seq.getSequence());
							hammerMap.put(genomeID, hammerSet);
						}
					}
				}
				log.info("{} genomes found with {}.", dnaMap.size(), sourName);
			}
			// Now we loop through the map, performing comparisons. This is a quadratic operation, since
			// we compare each instance to all the instances after it in the map (which is why the
			// NavigableMap interface is being used).
			log.info("Processing {} comparisons for {}.", dnaMap.size() * (dnaMap.size() - 1) / 2, sourName);
			this.genomeSet.stream().parallel().forEach(x -> this.runComparisons(sourName, x, dnaMap, hammerMap));
			log.info("Comparisons completed for {}.", sourName);
			this.reporter.finishSourReport();
		}
		log.info("Finishing report.");
		this.reporter.finishReport();
	}

	/**
	 * Process all the comparisons against the specified genome ID. If the genome ID
	 * has a map entry, its entry will be compared to every entry following it in the
	 * map. This is a concurrent method, so it synchronizes during reporting.
	 *
	 * @param sourName		name of the relevant SOUR
	 * @param hammerDb		hammer database
	 * @param genomeID		ID of the first genome to use in all the comparisons
	 * @param dnaMap		navigable map of genome IDs to DNA kmer objects
	 */
	private void runComparisons(String sourName, String genomeID, NavigableMap<String, DnaKmers> dnaMap,
			Map<String, Set<String>> hammerMap) {
		// Get the DNA kmers for the source genome.
		DnaKmers sourceKmers = dnaMap.get(genomeID);
		// Get the hammer set for the source genome.
		Set<String> sourceHammers = hammerMap.get(genomeID);
		// Only proceed if the genome had a sequence for this SOUR.
		if (sourceKmers != null) {
			// Loop through the entries for genomes whose ID is greater, saving comparison results.
			NavigableMap<String, DnaKmers> tailMap = dnaMap.tailMap(genomeID, false);
			List<HammerDnaDistReporter.Comparison> resultList = new ArrayList<>(tailMap.size());
			for (var tailEntry : tailMap.entrySet()) {
				String genome2 = tailEntry.getKey();
				// Compute the DNA distance.
				DnaKmers targetKmers = tailEntry.getValue();
				double dnaDist = sourceKmers.distance(targetKmers);
				// Compute the hammer distance.
				Set<String> targetHammers = hammerMap.get(genome2);
				double simCount = sourceHammers.stream().filter(x -> targetHammers.contains(x)).count();
				double hammerDist;
				if (simCount == 0.0)
					hammerDist = 1.0;
				else
					hammerDist = 1.0 - simCount / (sourceHammers.size() + targetHammers.size() - simCount);
				// Save the comparison result.
				resultList.add(new HammerDnaDistReporter.Comparison(genomeID, genome2, hammerDist, dnaDist));
			}
			// Now process the results. Note we process all the results at once to avoid the overhead of going
			// in and out of sync hundreds of times. Each process can run through all its comparisons uninterrupted
			// and only pause at the end when output is necessary.
			synchronized (this.reporter) {
				resultList.stream().forEach(x -> this.reporter.processComparison(x));
			}
		}
	}

	@Override
	public CompareFilter getFilter() {
		return this.filter;
	}

}
