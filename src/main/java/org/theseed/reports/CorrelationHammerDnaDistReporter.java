/**
 *
 */
package org.theseed.reports;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

/**
 * This report accumulates all the data for each SOUR and then outputs correlation coefficients.
 * The data is overwhelmed by 1/1 and 0/0 items, so these are removed before the correlation is
 * computed. The report will include a count of the "non-trivial" distance pairs, as well as
 * Kolmogorov-Smirnov coefficients for each of the correlated datasets.
 *
 * @author Bruce Parrello
 *
 */
public class CorrelationHammerDnaDistReporter extends HammerDnaDistReporter {

	// FIELDS
	/** stats object for the hammer distances */
	private DescriptiveStatistics hammerStats;
	/** stats object for the DNA distances */
	private DescriptiveStatistics dnaStats;
	/** computation helper for K-S test */
    private KolmogorovSmirnovTest ksTester;


	public CorrelationHammerDnaDistReporter(IParms processor) {
		super(processor);
		// Create a K-S test helper.
		this.ksTester = new KolmogorovSmirnovTest();
	}

	@Override
	protected void writeHeader() {
		this.println("SOUR\tpairs\tnon-trivial\tpearson\tkendall\tspearman\tK-S hammer\tK-S dna");
	}

	@Override
	protected void startSour(String sour) {
		// Set up the stats objects.
		this.hammerStats = new DescriptiveStatistics();
		this.dnaStats = new DescriptiveStatistics();
	}

	@Override
	protected void processUsefulComparison(Comparison comparison) {
		// Add the new values to the descriptor arrays.
		this.hammerStats.addValue(comparison.getHammerDist());
		this.dnaStats.addValue(comparison.getDnaDist());
	}

	@Override
	protected void summarizeSour(String sour, int totalCount, int usefulCount) {
		// The remaining numbers only exist if there are two or more nontrivial pairs.
		double pNum, kNum, sNum, ksHammer, ksDna;
		if (usefulCount >= 2) {
			PearsonsCorrelation pearson = new PearsonsCorrelation();
			pNum = pearson.correlation(this.hammerStats.getValues(), this.dnaStats.getValues());
			KendallsCorrelation kendall = new KendallsCorrelation();
			kNum = kendall.correlation(this.hammerStats.getValues(), this.dnaStats.getValues());
			SpearmansCorrelation spearman = new SpearmansCorrelation();
			sNum = spearman.correlation(this.hammerStats.getValues(), this.dnaStats.getValues());
			ksHammer = computeKS(this.hammerStats);
			ksDna = computeKS(this.dnaStats);
		} else {
			// Here we default to perfect correlation but not a normal distribution. In our case,
			// this latter is very bad, as it means our sample is not diverse enough.
			pNum = 1.0;
			kNum = 1.0;
			sNum = 1.0;
			ksHammer = 0.0;
			ksDna = 0.0;
		}
		// Write the output line.
		this.println(sour + "\t" + totalCount + "\t" + usefulCount + "\t" + pNum
				+ "\t" + kNum + "\t" + sNum + "\t" + ksHammer + "\t" + ksDna);
	}

	/**
	 * Compute the Kolmogorov-Smirnov value for a distance array.
	 *
	 * @param stats
	 *
	 * @return the Kolmogorov-Smirnov value; in our case a high value is good
	 */
	private double computeKS(DescriptiveStatistics stats) {
		// We create a normal distribution with the same mean and standard deviation.
		double mean = stats.getMean();
		double sdev = stats.getStandardDeviation();
		NormalDistribution normal = new NormalDistribution(mean, sdev);
		// Run the K-S test to determine the chance our dataset matches the distribution.
		double retVal = this.ksTester.kolmogorovSmirnovStatistic(normal, stats.getValues());
		return retVal;
	}

	@Override
	public void finishReport() {
	}

}
