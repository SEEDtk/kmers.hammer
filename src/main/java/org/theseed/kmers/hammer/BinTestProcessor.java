/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.binning.BinBuilder;
import org.theseed.binning.BlastTaxonomyComputer;
import org.theseed.binning.HammerBinningRule;
import org.theseed.binning.HammerBinningRule.IParms;
import org.theseed.binning.MappingBinRule;
import org.theseed.binning.TaxonomyComputer;
import org.theseed.counters.CountMap;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.utils.BaseReportProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command reads a hammer analysis file for a metagenomic contig file and then splits its FASTA
 * file into smaller files, one for each bin.  These smaller files can then be analyzed for taxonomy
 * and submitted to RAST.
 *
 * The hammer analysis file is tab-delimited with no headers.  The first column is a contig ID.  The
 * subsequent columns are in pairs, the first containing a reference genome ID and the second a count.
 * The counts will always be in descending order.  The first reference genome ID will always indicate
 * the bin into which the contig should be placed.  Based on various criteria, the others may cause
 * the placement to be suppressed.  Many contigs will have no pairs present.
 *
 * The hammer analysis file will be read from the standard input.  The positional parameters are
 * the name of the contigs file and the name of the output directory.  The standard output will contain
 * a report listing the output FASTA file names and the proposed taxonomic ID and name for each.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -i	input hammer analysis file (if not STDIN)
 * -o	output report file (if not STDOUT)
 *
 * --minCovg	minimum coverage for a contig to be kept (default 4)
 * --minLen		minimum length for a contig to be kept (default 400)
 * --seedFasta	FASTA file to use for determining taxonomy
 * --maxE		maximum acceptable E-level for taxonomy blast hits (default 1e-10)
 * --minMatch	mimumum fraction of the protein length that must match in a taxonomy blast hit (default 0.8)
 * --rule		rule to use for selecting the bin (default MAX)
 * --clear		erase the output directory before beginning
 * --minDiff	minimum hit difference for MAX binning rule
 * --taxTable	flat file to use for determining taxonomy when BLAST fails
 *
 * @author Bruce Parrello
 *
 */
public class BinTestProcessor extends BaseReportProcessor implements IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(BinTestProcessor.class);
    /** binning rule to use */
    private HammerBinningRule rule;
    /** map of bin IDs to domains */
    private Map<String, String> domainMap;
    /** map of bin IDs to backup taxonomic results */
    private Map<String, TaxonomyComputer.Result> backupTaxMap;
    /** taxonomy computer */
    private TaxonomyComputer taxComputer;

    // COMMAND-LINE OPTIONS

    /** SEED FASTA file for computing taxonomy */
    @Option(name = "--seedFasta", metaVar = "PhenTrnaSyntAlph.fa", usage = "DNA FASTA file for computing taxonomy",
            required = true)
    private File seedFastaFile;

    /** input file (if not STDIN) */
    @Option(name = "--input", aliases = { "-i" }, metaVar = "hammerHits.tsv", usage = "input file containing hammer hit data")
    private File inFile;

    /** minimum coverage for an acceptable contig */
    @Option(name = "--minCovg", metaVar = "5.0", usage = "minimum coverage for an acceptable contig")
    private double minCovg;

    /** minimum length for an acceptable contig */
    @Option(name = "--minLen", metaVar = "100", usage = "minimum length for an acceptable contig")
    private int minLen;

    /** maximum E-value for a taxonomy match */
    @Option(name = "--maxE", metaVar = "1e-10", usage = "maximum e-value for a taxonomy match")
    private double maxE;

    /** minimum match percent for a taxonomy match */
    @Option(name = "--minMatch", metaVar = "75.0", usage = "minimum length percentage for a taxonomy match")
    private double minMatchPercent;

    /** binning rule to use */
    @Option(name = "--rule", usage = "binning rule for accepting a contig into a particular bin")
    private HammerBinningRule.Type ruleType;

    /** TRUE to erase the output directory before starting */
    @Option(name = "--clear", usage = "if specified, the output directory will be cleared before starting")
    private boolean clearFlag;

    /** minimum hit count different for MAX binning rule */
    @Option(name = "--minDiff", metaVar = "1", usage = "the minimum hit count difference for the MAX binning rule")
    private int minDiff;

    /** taxonomy table file for backup taxonomy computation */
    @Option(name = "--taxTable", metaVar = "PhenTrnaSyntAlph.tbl", usage = "backup taxonomy result information file",
            required = true)
    private File taxTableFile;

    /** contigs file name */
    @Argument(index = 0, metaVar = "contigsIn.fna", usage = "file of contigs to bin", required = true)
    private File contigsInFile;

    /** output directory name */
    @Argument(index = 1, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    @Override
    protected void setReporterDefaults() {
        this.clearFlag = false;
        this.inFile = null;
        this.minCovg = 4.0;
        this.minLen = 400;
        this.maxE = 1e-10;
        this.minMatchPercent = 80.0;
        this.ruleType = HammerBinningRule.Type.MAX;
        this.minDiff = 8;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {

        if (this.inFile != null && ! this.inFile.canRead())
            throw new FileNotFoundException("Input file " + this.inFile + " is not found or unreadable.");
        if (! this.contigsInFile.canRead())
            throw new FileNotFoundException("Contig file " + this.contigsInFile + " is not found or unreadable.");
        if (this.minCovg <= 0.0)
            throw new FileNotFoundException("Minimum coverage must be at least 0.0.");
        if (this.minLen <= 0.0)
            throw new FileNotFoundException("Minimum length must be at least 0.");
        if (! this.seedFastaFile.canRead())
            throw new FileNotFoundException("SEED FASTA file " + this.seedFastaFile + " is not found or unreadable.");
        if (this.maxE < 0)
            throw new ParseFailureException("Maximum e-value cannot be negative.");
        if (this.minMatchPercent > 100.0 || this.minMatchPercent < 0)
            throw new ParseFailureException("Minimum match percent must be between 0 and 100.");
        if (! this.taxTableFile.canRead())
            throw new FileNotFoundException("Taxonomy table file " + this.taxTableFile + " is not found or unreadable.");
        // Create the binning rule.
        this.rule = this.ruleType.create(this);
        // Set up the output directory.
        if (! this.outDir.exists()) {
            // Here we must create the output directory.
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        } else if (! this.outDir.isDirectory())
            throw new IOException("Invalid output directory " + this.outDir);
        else {
            log.info("Output will be to {}.", this.outDir);
            if (this.clearFlag) {
                log.info("Erasing files in {}.", this.outDir);
                FileUtils.cleanDirectory(this.outDir);
            }
        }
        // Create the main taxonomy computer.
        log.info("Creating taxonomy computer from {}.", this.seedFastaFile);
        BlastParms parms = new BlastParms().maxE(this.maxE).pctLenOfSubject(this.minMatchPercent);
        this.taxComputer = new BlastTaxonomyComputer(this.seedFastaFile, parms);
        // Create the backup taxonomy tables.
        log.info("Getting backup taxonomy data from {}.", this.taxTableFile);
        this.backupTaxMap = new HashMap<String, TaxonomyComputer.Result>(8000);
        this.domainMap = new HashMap<String, String>(8000);
        try (TabbedLineReader taxStream = new TabbedLineReader(this.taxTableFile)) {
            int idCol = taxStream.findField("genome_id");
            int domCol = taxStream.findField("domain");
            int nameCol = taxStream.findField("rep_name");
            for (TabbedLineReader.Line line : taxStream) {
                String binId = line.get(idCol);
                String domain = line.get(domCol);
                this.domainMap.put(binId, domain);
                var result = new TaxonomyComputer.Result(StringUtils.substringBefore(binId, "."), line.get(nameCol));
                this.backupTaxMap.put(binId, result);
            }
            log.info("{} bin IDs set up for backup taxonomy.", this.backupTaxMap.size());
        }
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // First, we need to process the hammer input file to create the contig map.
        Map<String, String> contigMap = new HashMap<String, String>(1000);
        InputStream hammerStream;
        if (this.inFile == null) {
            hammerStream = System.in;
            log.info("Reading hammer data from standard input.");
        } else {
            hammerStream = new FileInputStream(this.inFile);
            log.info("Reading hammer data from {}.", this.inFile);
        }
        try (LineReader inStream = new LineReader(this.inFile)) {
            // Loop through the lines, computing bins.
            for (String line : inStream) {
                String[] cols = StringUtils.split(line, '\t');
                if (cols.length >= 3) {
                    // Here we have at least one hammer hit to process.
                    var counts = new CountMap<String>();
                    for (int i = 1; i < cols.length; i += 2)
                        counts.count(cols[i], Integer.valueOf(cols[i+1]));
                    // Compute the bin from the counts.
                    String binId = this.rule.choose(counts);
                    if (binId != null)
                        contigMap.put(cols[0], binId);
                }
            }
        } finally {
            // If we are reading hammer data from a file, close the file stream.
            if (this.inFile != null)
                hammerStream.close();
        }
        log.info("{} contigs had bins assigned.", contigMap.size());
        // Create the binning rule.
        var binRule = new MappingBinRule(contigMap);
        binRule.setMinCovg(this.minCovg);
        binRule.setMinLen(this.minLen);
        log.info("Minimum coverage is {} and minimum contig length is {}.", this.minCovg, this.minLen);
        // Set up the bin builder and open the contig input file.
        try (BinBuilder binner = new BinBuilder(binRule, this.outDir);
                FastaInputStream contigStream = new FastaInputStream(this.contigsInFile)) {
            // Read the contigs and make bin assignments.
            int contigsIn = 0;
            int contigsKept = 0;
            for (Sequence seq : contigStream) {
                contigsIn++;
                String bin = binner.binContig(seq);
                if (bin != null)
                    contigsKept++;
                if (log.isInfoEnabled() && contigsIn % 1000 == 0)
                    log.info("{} contigs processed, {} binned.", contigsIn, contigsKept);
            }
            log.info("{} contigs processed, {} binned.", contigsIn, contigsKept);
            // Close the binner to flush the output files.  Closing a closed binner is always safe,
            // so the autoclose will still work.
            binner.close();
            // Now we need to produce the output report about the bins.  This involves computing the
            // taxonomy of each bin, a nontrivial task.
            // Prime the output report.
            writer.println("bin_id\tfile\tcontigs\tdna\ttaxon\tdomain\tscientific_name");
            // Loop through the bins.
            var allBins = binner.getNonEmptyBins();
            for (BinBuilder.Stats bin : allBins) {
                // Get the bin FASTA file.
                File binFile = bin.getLocation();
                log.info("Processing bin FASTA file {}.", binFile);
                String binId = bin.getId();
                // Default to a virtual-bin case.
                TaxonomyComputer.Result taxResult = TaxonomyComputer.VIRTUAL_GENOME;
                String domain = "";
                if (! bin.isVirtual()) {
                    taxResult = taxComputer.analyzeFasta(binFile);
                    domain = this.domainMap.get(binId);
                    if (taxResult.equals(TaxonomyComputer.UNKNOWN_TAXON))
                        taxResult = this.backupTaxMap.getOrDefault(binId, TaxonomyComputer.UNKNOWN_TAXON);
                }
                // Compute the domain.
                writer.format("%s\t%s\t%d\t%d\t%s\t%s\t%s%n", binId, binFile.toString(), bin.getNum(),
                        bin.getLength(), taxResult.getId(), domain, taxResult.getName());
            }
        }

    }

    @Override
    public int getMinDiff() {
        return this.minDiff;
    }

}
