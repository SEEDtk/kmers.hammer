/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Spliterator;
import java.util.function.Consumer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;

/**
 * This object manages a directory containing sequences.  It allows iteration through the
 * members of the directory:  each member is either a fasta file or a genome file.  The
 * filename extension (.gto or .fna/.ffn/.fa/.fasta) determines the file type.  Each member
 * returned has an ID, a name, and a sequence list.  The object supports parallel streaming.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceDirectory implements Iterable<SequenceDirectory.Member> {

    /**
     * This is a late-binding parallelizing spliterator for processing the members in a stream.
     */
    public class Splitter implements Spliterator<Member> {

        /** array of files */
        private File[] files;
        /** current position of this splitter */
        private int pos;
        /** limit of this splitter */
        private int fence;

        /**
         * Construct a splitter for a portion of the file list.
         *
         * @param fileArray		array of files
         * @param start			starting index
         * @param end			limit index (past allocated end)
         */
        public Splitter(File[] fileArray, int start, int end) {
            this.files = fileArray;
            this.pos = start;
            this.fence = end;
        }

        @Override
        public boolean tryAdvance(Consumer<? super Member> action) {
            boolean retVal = false;
            // If there are more files in our section, return a member for the next one.
            if (this.pos < this.fence) {
                Member nextMember = SequenceDirectory.getMember(this.files[pos]);
                action.accept(nextMember);
                this.pos++;
                retVal = true;
            }
            return retVal;
        }

        @Override
        public Spliterator<Member> trySplit() {
            Spliterator<Member> retVal = null;
            // Attempt to split off a section from this spliterator.
            int splitPosition = (this.pos + this.fence) / 2;
            // Here we are too small to split.
            if (splitPosition <= this.pos)
                retVal = null;
            else {
                // Create a new splitter for the second section.
                retVal = new Splitter(this.files, splitPosition, this.fence);
                // Shrink us to the first section.
                this.fence = splitPosition;
            }
            return retVal;
        }

        @Override
        public long estimateSize() {
            // Return the number of items left.
            return this.fence - this.pos;
        }

        @Override
        public int characteristics() {
            return Spliterator.SUBSIZED + Spliterator.SIZED + Spliterator.IMMUTABLE + Spliterator.NONNULL;
        }

    }

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SequenceDirectory.class);
    /** array of sequence files in the directory */
    private final File[] files;
    /** pattern for extracting genome IDs from FASTA labels and comments */
    private static final Pattern GENOME_ID_FINDER = Pattern.compile("\\b(\\d+\\.\\d+)");
    /** pattern for parsing genome IDs from file names */
    private static final Pattern GENOME_ID_FILE = Pattern.compile("(\\d+\\.\\d+)\\..+");

    /**
     * This class acts as a file name filter for locating the sequence files in the input directory.
     */
    private static class SeqFileFilter implements FilenameFilter {

        @Override
        public boolean accept(File dir, String name) {
            String extension = StringUtils.substringAfterLast(name, ".");
            boolean retVal = false;
            switch (extension) {
            case "gto" :
            case "fna" :
            case "ffn" :
            case "fa" :
            case "fasta" :
                retVal = true;
            }
            return retVal;
        }

    }

    /**
     * This class represents a directory member returned from an iteration.
     */
    public static class Member {

        /** ID of the member (intended to be a genome ID) */
        private String id;
        /** name of the member (intended to be a genome or file name) */
        private String name;
        /** sequences in the file */
        private List<Sequence> sequences;
        /**
         * Construct a new directory member.
         *
         * @param genome_id		ID to use
         * @param genome_name	name to use
         * @param genome_seqs	sequences in the member
         */
        public Member(String genome_id, String genome_name, List<Sequence> genome_seqs) {
            this.id = genome_id;
            this.name = genome_name;
            this.sequences = genome_seqs;
        }
        /**
         * @return the member id
         */
        public String getId() {
            return this.id;
        }
        /**
         * @return the member name
         */
        public String getName() {
            return this.name;
        }
        /**
         * @return the sequence list
         */
        public List<Sequence> getSequences() {
            return this.sequences;
        }

    }

    /**
     * This class is the iterator for members of this sequence directory.
     */
    public class Iter implements Iterator<Member> {

        /** current position in file list */
        private int nextFileIdx;

        /**
         * Construct an iterator through the members of this sequence directory.
         */
        public Iter() {
            this.nextFileIdx = 0;
        }

        @Override
        public boolean hasNext() {
            return (this.nextFileIdx < SequenceDirectory.this.files.length);
        }

        @Override
        public Member next() {
            if (this.nextFileIdx >= SequenceDirectory.this.files.length)
                throw new NoSuchElementException("Attempt to advance past end of member file list.");
            // Get the next file.
            File memberFile = SequenceDirectory.this.files[this.nextFileIdx];
            Member retVal = SequenceDirectory.getMember(memberFile);
            // Increment for the next file.
            this.nextFileIdx++;
            // Return the member.
            return retVal;
        }


    }

    /**
     * Construct a sequence directory object for a given input directory.
     *
     * @param inDir		input directory to use for construction
     */
    public SequenceDirectory(File inDir) {
        // Save an array of the files containing sequence data.
        this.files = inDir.listFiles(new SeqFileFilter());
    }

    /**
     * @return the member corresponding to a file
     *
     * @param memberFile	file to convert into a member
     */
    protected static Member getMember(File memberFile) {
        Member retVal = null;
        try {
            // Determine the file type.
            String fileName = memberFile.getName();
            if (fileName.endsWith(".gto")) {
                // Here we have a genome file.
                Genome genome = new Genome(memberFile);
                retVal = new Member(genome.getId(), genome.getName(), genome.getSequences());
            } else {
                // Here we have a FASTA file.
                String id = "(unknown)";
                // Get the set of sequences from the file.
                List<Sequence> seqs = FastaInputStream.readAll(memberFile);
                // Try to find a genome ID in the header for the first sequence.
                if (seqs.size() >= 1) {
                    Sequence seq0 = seqs.get(0);
                    String firstLine = seq0.getLabel() + " " + seq0.getComment();
                    Matcher m = GENOME_ID_FINDER.matcher(firstLine);
                    if (m.find())
                        id = m.group(1);
                    else {
                        // Here we check for a file-name ID.
                        m = GENOME_ID_FILE.matcher(fileName);
                        if (m.matches())
                            id = m.group(1);
                    }
                }
                // Form the member.
                retVal = new Member(id, fileName, seqs);
            }
        } catch (IOException e) {
            // Percolate the IO exception as unchecked for compatability with iterators.
            throw new UncheckedIOException(e);
        }
        return retVal;
    }

    @Override
    public Iterator<Member> iterator() {
        return this.new Iter();
    }

    /**
     * @return the number of members in this directory
     */
    public int size() {
        return this.files.length;
    }

    /**
     * The spliterator is necessary to support parallel streams.  It is unordered, but sized.
     *
     * @return a spliterator for this directory
     */
    public Spliterator<Member> spliterator() {
        return new Splitter(this.files, 0, this.files.length);
    }

    /**
     * This method returns all the members as a parallel stream.
     */
    public Stream<Member> parallelStream() {
        return StreamSupport.stream(this.spliterator(), true);
    }

    /**
     * This method returns all the members as a serial stream.
     */
    public Stream<Member> stream() {
        return StreamSupport.stream(this.spliterator(), false);
    }

    /**
     * This method returns all the members as a stream with optional parallelization.
     * It provides an easy way to parameterize the parallelism.
     *
     * @param para		TRUE to parallelize the stream, else FALSE
     */
    public Stream<Member> stream(boolean para) {
        return StreamSupport.stream(this.spliterator(), para);
    }

}
