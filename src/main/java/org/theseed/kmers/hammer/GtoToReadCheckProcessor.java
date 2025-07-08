package org.theseed.kmers.hammer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3Genome;
import org.theseed.sequence.fastq.SampleDescriptor;
import org.theseed.sequence.fastq.SeqRead;
import org.theseed.sequence.GenomeKmers;
import org.theseed.sequence.fastq.FastqSampleGroup;
import org.theseed.sequence.fastq.ReadStream;
import org.theseed.utils.BasePipeProcessor;
import org.theseed.utils.FloatList;

/**
 * This command loads a GTO and memorizes its DNA kmers. Then it reads a FASTQ sample and finds
 * all the kmers in the sample that are not in the GTO. For each sample, we will report the
 * total number of kmers, the number that are not in the GTO, and how many good kmers and bad
 * kmers occur at different kmer quality levels. The quality of a kmer is expressed as the
 * change the kmer is correct. In addition, we will reject any read with an expected error of at least 1
 * nucleotide (though that can be changed by a command-line option).
 * 
 * The standard input should contain a tab-delimited file with a column containing genome IDs and a column
 * containing sample IDs. The positional parameters should be the name of a genome source from which the
 * GTOs will be loaded and the name of a directory containing the FASTQ samples. Each FASTQ sample must
 * be either a single file or a directory containing paired-end files. 
 * 
 * The output will be tab-delimited report on the standard outout.
 * 
 * The command-line options are as follows:
 * 
 * -h   display command-line usage
 * -v   display more frequent log messages
 * -i   input file name (default is standard input)
 * -o   output file name (default is standard output)
 * -K   kmer size (default 20)
 * 
 * --sampCol    index (1-based) or name of the column in the input file containing sample IDs (default "sample_id")
 * --idCol      index (1-based) or name of the column in the input file containing genome IDs (default "genome_id")
 * --source     type of genome source (default "DIR")
 * --filter     maximum number of expected bad base pairs for an acceptable read (default 1.0)
 * --quals      comma-delimited list of quality score cutoffs for reporting (default "0.7,0.8,0.9,0.95")
 * 
 * @author  Bruce Parrello
 */
public class GtoToReadCheckProcessor extends BasePipeProcessor {

    // FIELDS
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(GtoToReadCheckProcessor.class);
    /** genome source */
    private GenomeSource genomes;
    /** sample column index */
    private int sampColIdx;
    /** genome ID column index */
    private int gColIdx;
    /** FASTQ sample group */
    private FastqSampleGroup sampleGroup;
    /** array of hit quality cutoffs */
    private double[] qualCutoffs;

    // COMMAND-LINE OPTIONS

    /** index (1-based) or name of the column containing the sample IDs */
    @Option(name = "--sampCol", metaVar = "sample", usage = "index (1-based) or name of the column in the input file containing sample IDs")
    private String sampCol;

    /** index (1-based) or name of the column containing the genome IDs */
    @Option(name = "--idCol", metaVar = "genome", usage = "index (1-based) or name of the column in the input file containing genome IDs")
    private String gCol;

    /** type of genome source */
    @Option(name = "--source", usage = "type of genome source")
    private GenomeSource.Type sourceType;

    /** kmer size */
    @Option(name = "--kmerSize", aliases = { "-K" }, metaVar = "16", usage = "kmers size to use")
    private int kmerSize;

    /** maximum number of expected errors allowed in an acceptable read */
    @Option(name = "--filter", metaVar = "2.0", usage = "number of expected errors at or beyond which a read is rejected")
    private double badBaseFilter;

    /** list of quality score cutoffs for reporting */
    @Option(name = "--quals", metaVar = "0.2,0.3,0.5", usage = "comma-delimited list of quality score cutoffs for reporting")
    private String qualString;

    /** genome source name (file or directory) */
    @Argument(index = 0, metaVar = "genomeDir", usage = "name of the genome source (file or directory)", required = true)
    private File genomeDir;

    /** FASTQ sample directory */
    @Argument(index = 1, metaVar = "fastqDir", usage = "name of the directory containing FASTQ samples", required = true)
    private File fastqDir;

    @Override
    protected void setPipeDefaults() {
        this.sampCol = "sample_id";
        this.gCol = "genome_id";
        this.sourceType = GenomeSource.Type.DIR;
        this.kmerSize = 20;
        this.badBaseFilter = 1.0;
        this.qualString = "0.7,0.8,0.9,0.95";
    }

    @Override
    protected void validatePipeInput(TabbedLineReader inputStream) throws IOException {
        // Insure the input columns exist.
        this.sampColIdx = inputStream.findColumn(this.sampCol);
        this.gColIdx = inputStream.findColumn(this.gCol);
    }

    @Override
    protected void validatePipeParms() throws IOException, ParseFailureException {
        // Validate the kmer size.
        if (this.kmerSize < 2)
            throw new ParseFailureException("Kmer size must be at least 2.");
        // Validate the read quality filter.
        if (this.badBaseFilter < 0.0)
            throw new ParseFailureException("Bad base filter must be non-negative.");
        // Validate the quality cutoffs.
        FloatList myQualCutoffs = new FloatList(this.qualString);
        this.qualCutoffs = new double[myQualCutoffs.size() + 1];
        this.qualCutoffs[0] = 0.0;  // Always include 0.0 for the no-quality case.
        int i = 1;
        for (double q : myQualCutoffs) {
            if (q < 0.0 || q > 1.0)
                throw new ParseFailureException("Quality cutoff " + q + " is out of range [0,1].");
            this.qualCutoffs[i] = q;
            i++;
        }
        // Set the kmer size in the genome kmers object.
        GenomeKmers.setKmerSize(this.kmerSize);
        // Insure the genome source exists.
        if (! this.genomeDir.exists())
            throw new FileNotFoundException("Genome source directory " + this.genomeDir + " does not exist.");
        // Insure the FASTQ sample directory is valid.
        if (! this.fastqDir.isDirectory())
            throw new FileNotFoundException("FASTQ sample directory " + this.fastqDir + " is not found or invalid.");
        // Connect to the genome source.
        log.info("Connecting to genome source {} of type {}.", this.genomeDir, this.sourceType);
        this.genomes = this.sourceType.create(this.genomeDir);
        // Get the FASTQ sample group.
        log.info("Connecting to FASTQ sample directory {}.", this.fastqDir);
        this.sampleGroup = FastqSampleGroup.Type.FASTQ.create(this.fastqDir);
        log.info("{} samples in {}. {} genomes in {}.", this.sampleGroup.size(), this.fastqDir, this.genomes.size(), this.genomeDir);
    }

    @Override
    protected void runPipeline(TabbedLineReader inputStream, PrintWriter writer) throws Exception {
        // Create and write the report header.
        writer.print("sample_id\tgenome_id\tgenome_name\ttotal_reads\tbad_reads\tpct_good_reads\ttotal_kmers");
        for (double q : this.qualCutoffs) {
            String qString = String.format("%.2f", q);
            writer.print("\tq" + qString + "_good_kmers\tq" + qString + "_bad_kmers\tq" + qString + "_pct_good_kmers");
        }
        writer.println();
        // Set up counters for the number of lines read, the number of samples not found, and the number of samples processed.
        int inCount = 0;
        int procCount = 0;
        int skipCount = 0;
        // Loop through the input file, processing each sample.
        log.info("Reading sample/genome pairs from input.");
        for (TabbedLineReader.Line line : inputStream) {
            // Get the sample and genome IDs.
            String sampleId = line.get(this.sampColIdx);
            String genomeId = line.get(this.gColIdx);
            inCount++;
            // Get the sample.  If it is not found, we will skip this genome.
            SampleDescriptor sampleDesc = this.sampleGroup.getDescriptor(sampleId);
            if (sampleDesc == null) {
                log.warn("Sample {} not found in FASTQ directory {}.", sampleId, this.fastqDir);
                skipCount++;
            } else {
                // Get the genome.  If it is not found, we will skip this sample.
                Genome genome = this.genomes.getGenome(genomeId, P3Genome.Details.CONTIGS);
                if (genome == null) {
                    log.warn("Genome {} not found in source {}.", genomeId, this.genomeDir);
                    skipCount++;
                } else {
                    procCount++;
                    log.info("Processing sample {} {} against genome {}.", procCount, sampleId, genome);
                    long lastMsg = System.currentTimeMillis();
                    // Get the genome Kmers.
                    GenomeKmers genomeKmers = new GenomeKmers(genome);
                    // Now we read through the sample, counting the kmers. We need statistics objects for the
                    // quality scores. Note that we will not count unique kmers here, since the quality can vary
                    // between occurrences.
                    long kmerCount = 0;
                    // We also count the reads and how many pass the expected-error filter.
                    int readCount = 0;
                    int badReadCount = 0;
                    // Create counters for each quality cutoff.
                    long[] goodKmerCounts = new long[this.qualCutoffs.length];
                    long[] badKmerCounts = new long[this.qualCutoffs.length];
                    // Read through the sample.
                    try (ReadStream fastqStream = sampleDesc.reader()) {
                        for (SeqRead read : fastqStream) {
                            readCount++;
                            if (read.getExpectedErrors() >= this.badBaseFilter) {
                                // Skip this read: it has too many expected errors.
                                badReadCount++;
                            } else {
                                // Process the left sequence.
                                kmerCount += this.processReadPart(genomeKmers, read.getLseq(), read.getLQual(), goodKmerCounts, badKmerCounts);
                                // Process the right sequence.
                                kmerCount += this.processReadPart(genomeKmers, read.getRseq(), read.getRQual(), goodKmerCounts, badKmerCounts);
                                long now = System.currentTimeMillis();
                                if (now - lastMsg > 5000) {
                                    // Report progress every 10 seconds.
                                    log.info("Processed {} reads in sample {}.", readCount, sampleId);
                                    lastMsg = now;
                                }   
                            }    
                        }
                        // Write the results for this sample.
                        log.info("Finished processing sample {}: {} kmers processed in {} reads.", sampleId, kmerCount, readCount);
                        double pctGoodReads = (readCount == 0 ? 0.0 : (double) (readCount - badReadCount) * 100.0 / readCount);
                        writer.print(sampleId + "\t" + genomeId + "\t" + genome.getName() + "\t" + readCount + "\t" + badReadCount + "\t" 
                                + pctGoodReads + "\t" + kmerCount);
                        for (int i = 0; i < this.qualCutoffs.length; i++) {
                            long goodCount = goodKmerCounts[i];
                            long badCount = badKmerCounts[i];
                            double pctGood = (goodCount + badCount == 0 ? 0.0 : (double) goodCount * 100.0 / (goodCount + badCount));
                            writer.print("\t" + goodCount + "\t" + badCount + "\t" + pctGood);
                        }
                        writer.println();
                        writer.flush();
                    }
                }
            }   
        }
        log.info("Finished processing {} samples. {} lines read. {} samples skipped due to missing genomes or FASTQ files.", procCount, inCount, skipCount);
    }

    /**
     * This method processes a single read part (either left or right) and counts the good and bad kmers, plus
     * updates the statistics for the quality scores. A kmer is considered good if it is found in the genome kmers.
     * 
     * @param kmers         genome kmers to check against
     * @param seq           sequence to get the kmers from
     * @param qual          quality string for the sequence
     * @param goodStats     good-kmer counts for different quality score cutoffs
     * @param badStats      bad-kmer counts for different quality score cutoffs
     * 
     * @return the total number of kmers scanned in this read part
     */
    private int processReadPart(GenomeKmers kmers, String seq, String qual, long[] goodStats, long[] badStats) {
        // Loop through the kmers in the sequence.
        int retVal = 0;
        for (int i = 0; i <= seq.length() - this.kmerSize; i++) {
            // Get the kmer and its quality score.
            String kmer = seq.substring(i, i + this.kmerSize);
            double kQual = SeqRead.qualChance(qual, i, this.kmerSize);
            // Is the kmer good?
            long[] targetStats;
            if (kmers.contains(kmer)) {
                // Yes, it's good.
                targetStats = goodStats;
            } else {
                // No, it's bad.
                targetStats = badStats;
            }
            // Add the kmer to the total for the appropriate quality cutoffs.
            for (int j = 0; j < this.qualCutoffs.length; j++) {
                if (kQual >= this.qualCutoffs[j]) {
                    targetStats[j]++;
                }
            }
            retVal++;
        }
        return retVal;
    }

     
}
 
