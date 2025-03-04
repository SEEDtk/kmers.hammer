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
	public void openSourReport(String sour) {
		this.sourName = sour;
	}

	@Override
	public void processComparison(Comparison comparison) {
		this.println(this.sourName + "\t" + comparison.getGenome1() + "\t" + comparison.getGenome2()
				+ "\t" + comparison.getHammerDist() + "\t" + comparison.getDnaDist());
	}

	@Override
	public void finishSourReport() {
	}

	@Override
	public void finishReport() {
	}

}
