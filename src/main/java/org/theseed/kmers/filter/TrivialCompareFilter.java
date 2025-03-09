/**
 *
 */
package org.theseed.kmers.filter;

import org.theseed.reports.HammerDnaDistReporter.Comparison;

/**
 * This filter removes comparisons for which both distances are 1 or both are 0.
 *
 * @author Bruce Parrello
 *
 */
public class TrivialCompareFilter extends CompareFilter {

	public TrivialCompareFilter(IParms processor) {
		super(processor);
	}

	@Override
	public boolean isUseful(Comparison cmp) {
		double hammerDist = cmp.getHammerDist();
		double dnaDist = cmp.getDnaDist();
		return ! (hammerDist == dnaDist && (hammerDist == 1.0 || hammerDist == 0.0));
	}

}
