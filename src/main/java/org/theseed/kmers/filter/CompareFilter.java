/**
 *
 */
package org.theseed.kmers.filter;

import org.theseed.reports.HammerDnaDistReporter;

/**
 * This is the base class for comparison filters. Currently, it is used to eliminate
 * unusual comparisons that would otherwise distort the correlations. For example,
 * TRIVIAL removes comparisons where both distances are 1 or both are 0, while
 * EXTREME removes comparisons where either distance is 1 or 0.
 *
 * @author Bruce Parrello
 *
 */
public abstract class CompareFilter {

	/**
	 * This interface must be supported by the controlling command processor
	 * that uses this filter.
	 */
	public interface IParms {

	}

	/**
	 * This enum lists the types of filters supported.
	 */
	public static enum Type {
		/** remove comparisons where both distances are 1 or both are 0 */
		TRIVIAL {
			@Override
			public CompareFilter create(IParms processor) {
				return new TrivialCompareFilter(processor);
			}
		},
		/** remove comparisons where either distance is 1 or 0 */
		EXTREME {
			@Override
			public CompareFilter create(IParms processor) {
				return new ExtremeCompareFilter(processor);
			}
		};

		/**
		 * @return a comparison filter of this type
		 *
		 * @param processor		controlling command processor
		 */
		public abstract CompareFilter create(IParms processor);
	}

	/**
	 * Construct a new comparison filter.
	 *
	 * @param processor		controlling command processor
	 */
	public CompareFilter(IParms processor) {
	}

	/**
	 * @return TRUE if a comparison is useful, else FALSE
	 *
	 * @param cmp	comparison descriptor to check
	 */
	public abstract boolean isUseful(HammerDnaDistReporter.Comparison cmp);


}
