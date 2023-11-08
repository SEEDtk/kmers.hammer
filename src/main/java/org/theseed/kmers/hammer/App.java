package org.theseed.kmers.hammer;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * This program processes DNA hammers.  The commands are as follows.
 *
 *  hammers			find all the hammers for a specified set of genomes
 *  hammersF		use a protein-finder to find all the hammers for a specified repgen set
 *  closest			compute the closest genomes using the hammers
 *  contigs			create a FASTA file for a hammer processor from a genome directory
 *  define			create a hammer definition file from a roles.in.subsystems file
 *  dbLoad			load a hammer database from a flat file
 *  binTest			separate contigs into bins using hammers
 *  contigTest		analyze a contig FASTA file using hammers
 *  readTest		test a FASTQ sample and output unexpected hammer hits
 *  scanTest		analyze the results of a synthetic hammer test
 *  scanLocs		compare features hit by bad hammer matches in a synthetic test
 *  distTest		compute the ANI distances to the expected and actual choices in a synthetic test
 *  testStats		produce statistics for each test/match genome pair in a contig test
 *  sampReport		produce a bin report for a sample group using hammers
 *  gtoReport		count hammer hits in the genomes of a genome source
 *  pseudoBins		produce a simulated bin report from a binned sample group
 *  krakenBins		produce a simulated bin report from kraken output in a sample group
 *  synthBins		produce a simulated bin report from a synthetic contig file
 *  binComp			do a comparison of bin reports
 *  filter			filter a hammer database by removing hammers in other genomes
 *  report			produce a report on a hammer database
 *  hammerStats		determine which hammers were found in a contig-hammer run and how many of each
 *  debugMeta		isolate bad hits in a synthetic-sample contigTest used for abundance measures
 *  sourMap			create a map between hammer fids and SOUR role IDs
 *  hCompare		compare two hammer databases
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        BaseProcessor processor;
        switch (command) {
        case "hammers" :
            processor = new HammerProcessor();
            break;
        case "hammersF" :
            processor = new HammerFinderProcessor();
            break;
        case "closest" :
            processor = new FindClosestProcessor();
            break;
        case "define" :
            processor = new HammerDefineProcessor();
            break;
        case "contigs" :
            processor = new HammerFastaProcessor();
            break;
        case "dbLoad" :
            processor = new HammerDbLoadProcessor();
            break;
        case "binTest" :
            processor = new BinTestProcessor();
            break;
        case "contigTest" :
            processor = new ContigTestProcessor();
            break;
        case "readTest" :
            processor = new ReadTestProcessor();
            break;
        case "scanTest" :
            processor = new ContigTestAnalysisProcessor();
            break;
        case "scanLocs" :
            processor = new ContigTestLocationProcessor();
            break;
        case "distTest" :
            processor = new ContigTestDistanceProcessor();
            break;
        case "testStats" :
            processor = new ContigTestStatisticsProcessor();
            break;
        case "sampReport" :
            processor = new SampleBinReportProcessor();
            break;
        case "pseudoBins" :
            processor = new PseudoBinReportProcessor();
            break;
        case "krakenBins" :
            processor = new KrakenBinReportProcessor();
            break;
        case "synthBins" :
            processor = new SynthBinReportProcessor();
            break;
        case "binComp" :
            processor = new BinReportCompareProcessor();
            break;
        case "filter" :
            processor = new HammerFilterProcessor();
            break;
        case "report" :
            processor = new HammerReportProcessor();
            break;
        case "hammerStats" :
            processor = new HammerStatsProcessor();
            break;
        case "debugMeta" :
            processor = new DebugMetaProcessor();
            break;
        case "sourMap" :
            processor = new HammerSourMapProcessor();
            break;
        case "gtoReport" :
            processor = new GtoHammerReportProcessor();
            break;
        case "hCompare" :
            processor = new HammerCompareReportProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
