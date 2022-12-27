package org.theseed.kmers.hammer;

import java.util.Arrays;

import org.theseed.utils.BaseProcessor;

/**
 * This program processes DNA hammers.  The commands are as follows.
 *
 *  hammers		find all the hammers for a specified set of genomes
 *  hammersF	use a protein-finder to find all the hammers for a specified repgen set
 *  closest		compute the closest genomes using the hammers
 *  contigs		create a FASTA file for a hammer processor from a genome directory
 *  define		create a hammer definition file from a roles.in.subsystems file
 *  dbLoad		load a hammer database from a flat file
 *  binTest		separate contigs into bins using hammers
 *  contigTest	analyze a contig FASTA file using hammers
 *  scanTest	analyze the results of a synthetic hammer test
 *  scanLocs	compare features hit by bad hammer matches in a synthetic test
 *  distTest	compute the ANI distances to the expected and actual choices in a synthetic test
 *  testStats	produce statistics for each test/match genome pair in a contig test
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
        default :
            throw new RuntimeException("Invalid command " + command + ".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
