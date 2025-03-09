/**
 *
 */
package org.theseed.kmers.filter;

import org.theseed.reports.HammerDnaDistReporter.Comparison;

/**
 * This filter removes comparisons where either distance is 1 or 0.
 *
 * @author Bruce Parrello
 *
 */
public class ExtremeCompareFilter extends CompareFilter {

	public ExtremeCompareFilter(IParms processor) {
		super(processor);
	}

	@Override
	public boolean isUseful(Comparison cmp) {
		double hammerDist = cmp.getHammerDist();
		double dnaDist = cmp.getDnaDist();
		return ! (hammerDist == 1.0 || hammerDist == 0.0 || dnaDist == 1.0 || dnaDist == 0.0);
	}

}
