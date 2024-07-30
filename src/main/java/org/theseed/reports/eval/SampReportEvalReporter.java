/**
 *
 */
package org.theseed.reports.eval;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Set;

import org.theseed.counters.CountMap;
import org.theseed.reports.BaseWritingReporter;

/**
 * This is the base clase for all sampReport evaluation reports.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SampReportEvalReporter extends BaseWritingReporter {

    // FIELDS
    /** controlling command processor */
    private IParms processor;

    /**
     * This interface is used to insure the controlling processor has the parameters needed
     * to process the reports.
     */
    public interface IParms {

        /**
         * @return a descriptor for the specified sample
         *
         * @param sampleId	ID of the sample whose descriptor is desired
         */
        public SampleDescriptor getDescriptor(String sampleId);

        /**
         * @return the distance between two genomes, or NaN if unknown
         */
        public double getDistance(String genome1, String genome2);

        /**
         * @return the set of roles of interest
         */
        public Set<String> getRoleSet();

    }

    /**
     * This enumeration describes the implemented report types.
     */
    public static enum Type {
        /** listing of the best reports for each sample */
        RATINGS {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new RatingSampReportEvalReporter(processor, writer);
            }
        },
        /** summary of good and bad hits in each report */
        QUALITY {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new QualitySampReportEvalReporter(processor, writer);
            }
        },
        /** details of results for each report */
        DETAILS {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new DetailSampReportEvalReporter(processor, writer);
            }
        },
        /** compare good and bad hits per sample by role */
        ROLES {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new RoleSummarySampReportEvalReporter(processor, writer);
            }
        },
        /** summary of individual sample performance in each report */
        SAMPLE {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new SamplingSampReportEvalReporter(processor, writer);
            }
        },
        /** comparison of distances to hit percentage for each report */
        COMPARE {
            @Override
            public SampReportEvalReporter create(IParms processor, PrintWriter writer) {
                return new CompareSampReportEvalReporter(processor, writer);
            }
        };

        /**
         * Create a reporter of this type
         *
         * @param processor		controlling command processor
         * @param writer		output print writer
         *
         * @return a reporter of this type
         */
        public abstract SampReportEvalReporter create(IParms processor, PrintWriter writer);
    }

    /**
     * This class describes the data we need for a sample from the input file.
     */
    public static class SampleDescriptor {

        /** sample ID */
        private String sampleId;
        /** sample genome ID */
        private String genomeId;
        /** sample genome name */
        private String genomeName;
        /** expected representative ID */
        private String repId;

        /**
         * @return the sampleId
         */
        public String getSampleId() {
            return this.sampleId;
        }

        /**
         * @return the genomeId
         */
        public String getGenomeId() {
            return this.genomeId;
        }

        /**
         * @return the genomeName
         */
        public String getGenomeName() {
            return this.genomeName;
        }

        /**
         * @return the repId
         */
        public String getRepId() {
            return this.repId;
        }

        /**
         * Specify a new sample ID.
         *
         * @param sampleId 	the sampleId to set
         */
        public void setSampleId(String sampleId) {
            this.sampleId = sampleId;
        }

        /**
         * Specify a new genome ID.
         *
         * @param genomeId 	the genomeId to set
         */
        public void setGenomeId(String genomeId) {
            this.genomeId = genomeId;
        }

        /**
         * Specify a new genome name.
         *
         * @param genomeName 	the genomeName to set
         */
        public void setGenomeName(String genomeName) {
            this.genomeName = genomeName;
        }

        /**
         * Specify a new representative genome ID.
         *
         * @param repId 	the repId to set
         */
        public void setRepId(String repId) {
            this.repId = repId;
        }

    }

    /**
     * Create a reporter with a specified controlling command processor and a specified output writer.
     *
     * @param processor		controlling command processor
     * @param writer		output report writer
     */
    public SampReportEvalReporter(IParms processor, PrintWriter writer) {
        super(writer);
        this.processor = processor;
    }

    /**
     * Begin processing a new bin report.
     *
     * @param reportName	name of bin report
     */
    public abstract void openFile(String name);

    /**
     * Record hits against a sample.
     *
     * @param desc			relevant sample descriptor
     * @param repId			ID of the repgen hit
     * @param repName		name of the repgen hit
     * @param count			weighted hit count
     * @param roleCount 	number of roles hit
     * @param roleCounts 	map of number of hits per each role
     *
     * @throws IOException
     */
    public abstract void recordHits(SampleDescriptor desc, String repId, String repName, double count, int roleCount, CountMap<String> roleCounts) throws IOException;

    /**
     * Finish processing a bin report.
     *
     * @throws IOException
     */
    public abstract void closeFile() throws IOException;

    /**
     * Finish the entire report.
     */
    public abstract void closeReport();

    /**
     * @return the descriptor for a specified sample
     *
     * @param sampleId	ID of the relevant sample
     */
    public SampleDescriptor getSample(String sampleId) {
        return this.processor.getDescriptor(sampleId);
    }

    /**
     * @return the distance between two genomes, or NaN if it is unknown
     *
     * @param g1	ID of the first genome
     * @param g2	ID of the second genome
     */
    public double getDistance(String g1, String g2) {
        return this.processor.getDistance(g1, g2);
    }

    /**
     * @return the distance between a sample and a genome, or NaN if it is unknown
     *
     * @param sample	sample ID
     * @param genome	genome ID
     */
    public double getSampleDistance(SampleDescriptor sample, String genome) {
        String g1 = sample.getGenomeId();
        return this.processor.getDistance(g1, genome);
    }

}
