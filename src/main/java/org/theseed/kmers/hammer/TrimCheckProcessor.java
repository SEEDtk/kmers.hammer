/**
 *
 */
package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.fastq.FastqSampleDescriptor;
import org.theseed.sequence.fastq.ReadStream;
import org.theseed.sequence.fastq.SampleDescriptor;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.utils.BaseProcessor;
import org.theseed.utils.ParseFailureException;

/**
 * This command produces a FASTA file from the differences between trimmed and untrimmed FASTQ files.  The
 * trimmed files should have an extension of "_ptrim.fq" or "_ptrim.fastq" and the untrimmed files an extension
 * of ".fq" or ".fastq".
 *
 * The positional parameter should be the name of the input directory and the name of the sample.  The FASTQ file
 * names should all begin with the sample name.  If the file name contains "_2" or "_R2" it will be treated as a
 * reverse-end file; otherwise it will be treated as a forward-end file.
 *
 * The FASTA file will be produced on the standard output.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more frequent log information
 * -o	output FASTA file (if not STDOUT)
 *
 * --qual	minimum acceptable sequence quality (default 10.0)
 * --kmer	kmer size for hammers (default 20)
 *
 * @author Bruce Parrello
 *
 */
public class TrimCheckProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TrimCheckProcessor.class);
    /** trimmed sample */
    private SampleDescriptor trimmedSample;
    /** untrimmed sample */
    private SampleDescriptor untrimmedSample;
    /** output stream */
    private FastaOutputStream outStream;
    /** number of downgraded reads */
    private int downgrades;
    /** number of rejected reads */
    private int rejected;
    /** total quality of rejected reads */
    private double rejectedQual;
    /** maximum quality of a rejected read */
    private double maxRejectQual;
    /** number of adapters found */
    private int adaptersFound;
    /** number of untrimmed reads found */
    private int untrimCount;
    /** number of shortened sequences */
    private int shortened;
    /** number of bytes removed by shortening */
    private int shortAmount;
    /** number of sequences output */
    private int seqCount;
    /** file filter to return only possible FASTQ files */
    private static final FilenameFilter FASTQ_FILENAME_FILTER = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            boolean retVal = false;
            File pathname = new File(dir, name);
            if (pathname.isFile())
                retVal = (name.endsWith(".fq") || name.endsWith(".fastq"));
            return retVal;
        }

    };

    // COMMAND-LINE OPTIONS

    /** minimum acceptable quality for an input sequence */
    @Option(name = "--qual", metaVar = "30.0", usage = "minimum acceptable sequence quality (0 to 99)")
    private double minQual;

    /** kmer size for the hammers */
    @Option(name = "--kmer", metaVar = "21", usage = "kmer size used for hammers")
    private int kmerSize;

    /** output FASTA file name */
    @Option(name = "--output", aliases = { "-o" }, metaVar = "outFile.fna", usage = "output FASTA file name (if not STDOUT)")
    private File outFile;

    /** input directory name */
    @Argument(index = 0, metaVar = "inDir", usage = "input directory name", required = true)
    private File inDir;

    /** input sample name */
    @Argument(index = 1, metaVar = "sampleName", usage = "sample file name prefix", required = true)
    private String sampleName;

    @Override
    protected void setDefaults() {
        this.minQual = 10.0;
        this.kmerSize = 20;
        this.outFile = null;
        this.outStream = null;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Validate the kmer size.
        if (this.kmerSize < 1)
            throw new ParseFailureException("Kmer size must be at least 1.");
        // Verify that the input directory exists.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Look for the input files in the input directory.
        String[] possibleInFiles = this.inDir.list(FASTQ_FILENAME_FILTER);
        // We will save the file names in here.
        String[] trim = new String[] { null, null };
        String[] untrim = new String[] { null, null };
        for (String possibleInFile : possibleInFiles) {
            if (possibleInFile.startsWith(this.sampleName)) {
                // Here the file belongs to our sample.  Figure out which of the two samples it is.
                String[] targetPair;
                String type;
                if (possibleInFile.contains("_ptrim.f")) {
                    targetPair = trim;
                    type = "trimmed";
                } else {
                    targetPair = untrim;
                    type = "untrimmed";
                }
                // Now determine which direction it is.
                int idx = (possibleInFile.contains("_R2") || possibleInFile.contains("_2")) ? 1 : 0;
                targetPair[idx] = possibleInFile;
                log.info("Using {} as file {} for {} sample.", possibleInFile, (idx + 1), type);
            }
        }
        // Verify we have a file for each mode (trimmed, untrimmed).
        if (trim[0] == null && trim[1] == null || untrim[0] == null && untrim[1] == null)
            throw new ParseFailureException("Input directory does not contain both trimmed and untrimmed samples.");
        // Create the sample descriptors.
        log.info("Connecting to sample streams.");
        this.trimmedSample = new FastqSampleDescriptor(this.inDir, "trimmed", trim[0], trim[1]);
        this.untrimmedSample = new FastqSampleDescriptor(this.inDir, "raw", untrim[0], untrim[1]);
        // Now set up the output.
        if (this.outFile == null) {
            log.info("FASTA will be produced on the standard output.");
            this.outStream = new FastaOutputStream(System.out);
        } else {
            log.info("FASTA will be written to {}.", this.outFile);
            this.outStream = new FastaOutputStream(this.outFile);
        }
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // Initialize the counters.
        this.adaptersFound = 0;
        this.downgrades = 0;
        this.rejected = 0;
        this.seqCount = 0;
        this.untrimCount = 0;
        this.shortened = 0;
        this.shortAmount = 0;
        // Initialize the reject-quality trackers.
        this.maxRejectQual = 0.0;
        this.rejectedQual = 0.0;
        // Connect to the two read streams.
        try (ReadStream untrimmedReads = this.untrimmedSample.reader();
                ReadStream trimmedReads = this.trimmedSample.reader()) {
            // We know that every trimmed read occurs in the untrimmed sample, and that they are in the same order.
            // This allows us to read both files in parallel.  We keep the current trimmed read in memory, and loop
            // through the untrimmed reads.
            SeqRead trimmedCurr = trimmedReads.safeNext();
            for (SeqRead untrimmedCurr : untrimmedReads) {
                this.untrimCount++;
                String untrimmedId = untrimmedCurr.getLabel();
                double untrimmedQual = untrimmedCurr.getQual();
                if (this.checkReads(untrimmedId, trimmedCurr)) {
                    // Here we have a match for the trimmed read.  Get its information.
                    double trimmedQual = trimmedCurr.getQual();
                    // The sequence will appear in both streams.  First, check if it will be
                    // absent from the trimmed scan but not the untrimmed scan.  This will
                    // happen if the trimmer decides to downgrade the quality.
                    if (trimmedQual < this.minQual && untrimmedQual >= this.minQual) {
                        this.writeSequence(untrimmedId, untrimmedCurr.getSequence(), "downgraded", trimmedQual);
                        this.downgrades++;
                    } else {
                        // Check for an adapter trim.
                        this.checkTrim(untrimmedId, trimmedCurr.getLseq(), untrimmedCurr.getLseq(), "left");
                        this.checkTrim(untrimmedId, trimmedCurr.getRseq(), untrimmedCurr.getRseq(), "right");
                    }
                    // Push forward to the next trimmed read.
                    trimmedCurr = trimmedReads.safeNext();
                } else {
                    // Here the untrimmed sequence was rejected by the trimmer.
                    this.writeSequence(untrimmedId, untrimmedCurr.getSequence(), "rejected", untrimmedQual);
                    this.rejected++;
                    this.rejectedQual += untrimmedQual;
                    if (this.maxRejectQual < untrimmedQual)
                        this.maxRejectQual = untrimmedQual;
                }
                if (log.isInfoEnabled() && this.untrimCount % 50000 == 0)
                    log.info("{} untrimmed reads processed.  {} sequences output.", this.untrimCount, this.seqCount);
            }
            // Flush the output.
            this.outStream.flush();
            // Now we need to log the counters.
            log.info("{} untrimmed reads processed.", this.untrimCount);
            log.info("{} sequences output.  {} reads downgraded, {} rejected, {} adapters found.", this.seqCount, this.downgrades,
                    this.rejected, this.adaptersFound);
            // Do the metrics on the rejections.
            if (this.rejected > 0) {
                double meanReject = this.rejectedQual / this.rejected;
                log.info("Mean quality of rejected sequences was {}. Max was {}.", meanReject, this.maxRejectQual);
            }
            // Do the metrics on the shortening.
            if (this.shortened > 0) {
                double meanShort = this.shortAmount / (double) this.shortened;
                log.info("{} sequences shortened, mean amount {}.", this.shortened, meanShort);
            }
        } finally {
            // Insure the output file is closed.
            if (this.outFile != null)
                this.outStream.close();
        }

    }

    /**
     * @return TRUE if the untrimmed-read ID matches the ID of the specified trimmed read (which may be NULL)
     *
     * @param untrimmedId	ID of the current untrimmed read
     * @param trimmedCurr	current trimmed read (or NULL if there is none)
     */
    private boolean checkReads(String untrimmedId, SeqRead trimmedCurr) {
        boolean retVal = (trimmedCurr != null && untrimmedId.contentEquals(trimmedCurr.getLabel()));
        return retVal;
    }

    /**
     * Write a sequence to the output FASTA stream.
     *
     * @param untrimmedId	ID of the source untrimmed read
     * @param sequence		sequence to write
     * @param comment		comment describing the reason for the sequence
     * @param qual			quality of the sequence
     *
     * @throws IOException
     */
    private void writeSequence(String untrimmedId, String sequence, String comment, double qual) throws IOException {
        this.seqCount++;
        // Format the label.  We use a counter to insure it's unique.
        String label = String.format("%s.%08d", this.sampleName, this.seqCount);
        // Format the full comment.
        StringBuilder outComment = new StringBuilder(80);
        outComment.append(untrimmedId).append('\t').append(comment).append('\t');
        if (Double.isFinite(qual))
            outComment.append(qual);
        // Write the sequence.
        this.outStream.write(new Sequence(label, outComment.toString(), sequence));
    }

    /**
     * Check for an adapter removed from a left or right read sequence.
     *
     * @param untrimmedId		ID of the sequence
     * @param trimmedSeq		trimmed version of sequence
     * @param untrimmedSeq		untrimmed version of sequence
     * @param comment			comment to use during output
     *
     * @throws IOException
     */
    private void checkTrim(String untrimmedId, String trimmedSeq, String untrimmedSeq, String comment) throws IOException {
        // Check for a shortened sequence.  Note that a sequence may be NULL if it's the missing side of a singleton.
        int originalLen = (untrimmedSeq == null ? 0 : untrimmedSeq.length());
        int newLen = (trimmedSeq == null ? 0 : trimmedSeq.length());
        if (newLen < originalLen) {
            // Here base pairs were trimmed.
            this.shortened++;
            int diff = originalLen - newLen;
            this.shortAmount += diff;
            // This method returns FALSE if either sequence is NULL (which would be the case if the trimmed sequence
            // was completely deleted).
            if (StringUtils.startsWith(untrimmedSeq, trimmedSeq)) {
                // Find the adapter.  We enforce a minimum equal to the hammer size.
                int pos;
                if (diff < this.kmerSize)
                    pos = originalLen - this.kmerSize;
                else
                    pos = originalLen - diff;
                String adapter = untrimmedSeq.substring(pos);
                this.adaptersFound++;
                // Write the sequence.  We don't care about quality here.
                String outComment = String.format("%s trim length %d", comment, diff);
                this.writeSequence(untrimmedId, adapter, outComment, Double.NaN);
            }
        }
    }


}
