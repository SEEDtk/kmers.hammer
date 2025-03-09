/**
 *
 */
package org.theseed.reports;

/**
 * This report echos the raw distance data to the output.
 *
 * @author Bruce Parrello
 *
 */
public class ListHammerDnaDistReporter extends HammerDnaDistReporter {

	// FIELDS
	/** name of the current SOUR */
	private String sourName;

	public ListHammerDnaDistReporter(IParms processor) {
		super(processor);
	}

	@Override
	protected void writeHeader() {
		this.println("sour\tgenome1\tgenome2\thammers\tDNA");
	}

	@Override
	protected void startSour(String sour) {
		this.sourName = sour;
	}

	@Override
	protected void processUsefulComparison(Comparison comparison) {
		this.println(this.sourName + "\t" + comparison.getGenome1() + "\t" + comparison.getGenome2()
				+ "\t" + comparison.getHammerDist() + "\t" + comparison.getDnaDist());
	}

	@Override
	public void summarizeSour(String sour, int totalCount, int usefulCount) {
	}

	@Override
	public void finishReport() {
	}

}
