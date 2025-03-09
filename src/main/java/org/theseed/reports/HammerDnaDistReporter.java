/**
 *
 */
package org.theseed.reports;

import org.theseed.kmers.filter.CompareFilter;

/**
 * This is the base class for all reports comparing hammer and DNA distances.
 *
 * @author Bruce Parrello
 *
 */
public abstract class HammerDnaDistReporter extends BaseReporterReporter {

	// FIELDS
	/** comparison filter */
	private CompareFilter filter;
	/** total number of comparisons */
	private int total;
	/** number of useful comparisons */
	private int useful;
	/** name of current SOUR */
	private String sourName;

	/**
	 * This interface must be supported by any controlling command processor that
	 * generates a hammer/DNA distance comparison report.
	 */
	public interface IParms {

		/**
		 * @return the comparison filter
		 */
		public CompareFilter getFilter();

	}

	/**
	 * type of report
	 */
	public static enum Type {
		/** show correlation between hammer distance and DNA distance for each SOUR */
		CORRELATION {
			@Override
			public HammerDnaDistReporter create(IParms processor) {
				return new CorrelationHammerDnaDistReporter(processor);
			}
		},
		/** list the comparison results */
		LIST {
			@Override
			public HammerDnaDistReporter create(IParms processor) {
				return new ListHammerDnaDistReporter(processor);
			}
		};

		/**
		 * @return a reporter of this type
		 *
		 * @param processor		controlling command processor
		 */
		public abstract HammerDnaDistReporter create(IParms processor);

	}

	/**
	 * This is a useful utility class that contains the results of a comparison.
	 */
	public static class Comparison {

		/** ID of the first genome */
		private String genome1;
		/** ID of the second genome */
		private String genome2;
		/** hammer distance */
		private double hammerDist;
		/** DNA distance */
		private double dnaDist;

		/**
		 * Construct a comparison object.
		 *
		 * @param g1		ID of the first genome
		 * @param g2		ID of the second genome
		 * @param hammer	hammer distance
		 * @param dna		DNA distance
		 */
		public Comparison(String g1, String g2, double hammer, double dna) {
			this.genome1 = g1;
			this.genome2 = g2;
			this.hammerDist = hammer;
			this.dnaDist = dna;
		}

		/**
		 * @return the ID of the first genome
		 */
		public String getGenome1() {
			return this.genome1;
		}

		/**
		 * @return the ID of the second genome
		 */
		public String getGenome2() {
			return this.genome2;
		}

		/**
		 * @return the hammer distance
		 */
		public double getHammerDist() {
			return this.hammerDist;
		}

		/**
		 * @return the dna kmer distance
		 */
		public double getDnaDist() {
			return this.dnaDist;
		}

	}

	/**
	 * Construct a hammer/DNA distance comparison report.
	 *
	 * @param processor		controlling command processor
	 */
	public HammerDnaDistReporter(IParms processor) {
		this.filter = processor.getFilter();
	}

	/**
	 * Initialize the data structures for a new SOUR.
	 *
	 * @param sourName	name of the SOUR
	 */
	final public void openSourReport(String sour) {
		// Clear the counters.
		this.total = 0;
		this.useful = 0;
		// Save the SOUR name.
		this.sourName = sour;
		// Give the subclass a chance to initialize.
		this.startSour(sour);
	}

	/**
	 * Initialize the report section for a SOUR.
	 *
	 * @param sour	name of the SOUR
	 */
	protected abstract void startSour(String sour);

	/**
	 * Record the result for a single comparison.
	 *
	 * @param comparison	comparison result to record
	 */
	final public void processComparison(Comparison comparison) {
		this.total++;
		if (this.filter.isUseful(comparison)) {
			this.useful++;
			this.processUsefulComparison(comparison);
		}

	}

	/**
	 * Process a comparison that goes into the report.
	 *
	 * @param comparison	comparison to process
	 */
	protected abstract void processUsefulComparison(Comparison comparison);

	/**
	 * Total and optionally output data for a single SOUR.
	 */
	final public void finishSourReport() {
		this.summarizeSour(this.sourName, this.total, this.useful);
	}

	/**
	 * Produce the summary information for the current SOUR.
	 *
	 * @param sour			name of the SOUR
	 * @param totalCount	total number of comparisons
	 * @param usefulCount	number of useful comparisons
	 */
	protected abstract void summarizeSour(String sour, int totalCount, int usefulCount);

	/**
	 * Total and optionally output data for the entire report.
	 */
	public abstract void finishReport();

}
