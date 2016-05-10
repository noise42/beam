package BEAM2;

public class SummarizedData {
//Contiene metodi e campi per scrivere il summary di fine run
	/*Score
	MCC - SearchObj.computeMCC (return)
	score - multirun
	consensus - SearchObj.computeConsensus (return)
	coverage - TP/P (SearchObj)
	input file - contenuto in CLP
	input sequences size - P (SearchObj)
	background file - contenuto in CLP
	bg size - N (SearchObj)
	mask - lo si prende da MultiRun
	mean structure content (msc) - da calcolare per il bin (da dati iniziali!)
	mean length (ml) - da calcolare per il bin (da dati iniziali!)
	*/
	
	/*
	 * in particolare 
	 * 		inputFile, bgFile possono essere riempiti prima
	 * 
	 */
	private double MCC, coverage,fallout,acc, msc, ml, score;
	private int seqSize,bgSize, mask;
	private String consensus, qbearConsensus, inputFile, bgFile;
	
	SummarizedData(){
		MCC=coverage=msc=ml=score=fallout=acc=0.0;
		seqSize=bgSize=mask=0;
		consensus=inputFile=bgFile="";
	}
	//GETTERS
	
	public double getMCC() {
		return MCC;
	}
	public double getCoverage() {
		return coverage;
	}
	public double getFallout() {
		return fallout;
	}
	public double getAcc() {
		return acc;
	}
	public double getMsc() {
		return msc;
	}
	public double getMl() {
		return ml;
	}
	public double getScore() {
		return score;
	}
	public int getInputSize() {
		return seqSize;
	}
	public int getBgSize() {
		return bgSize;
	}
	public int getMask() {
		return mask;
	}
	public String getConsensus() {
		return consensus;
	}
	public String getqBEARConsensus() {
		return qbearConsensus;
	}
	public String getInputFile() {
		return inputFile;
	}
	public String getBgFile() {
		return bgFile;
	}
	
	//SETTERS
	
	public void setMCC(double mCC) {
		MCC = mCC;
	}
	public void setCoverage(double coverage) {
		this.coverage = coverage;
	}
	public void setFallout(double fallout) {
		this.fallout = fallout;
	}
	public void setAcc(double acc) {
		this.acc = acc;
	}
	public void setMsc(double msc) {
		this.msc = msc;
	}
	public void setMl(double ml) {
		this.ml = ml;
	}
	public void setScore(double score) {
		this.score = score;
	}
	public void setInputSize(int seqSize) {
		this.seqSize = seqSize;
	}
	public void setBgSize(int bgSize) {
		this.bgSize = bgSize;
	}
	public void setMask(int mask) {
		this.mask = mask;
	}
	public void setConsensus(String consensus) {
		this.consensus = consensus;
	}
	public void setqbearConsensus(String qconsensus) {
		this.qbearConsensus = qconsensus;
	}
	public void setInputFile(String inputFile) {
		this.inputFile = inputFile;
	}
	public void setBgFile(String bgFile) {
		this.bgFile = bgFile;
	}

	public String computeBin() {
		// restituisce LSbin di appartenenza
		String l="";
		String s="";
		if (this.ml < 88) l="1";
		else if(this.ml < 128) l="2";
		else if(this.ml <= 500) l="3";
		else l="4";

		if (this.msc < .48) s="1";
		else if(this.msc < .5555556) s="2";
		else if(this.msc < .627907) s="3";
		else if(this.msc < .7027027 ) s="4";
		else s="5";
	
		return l+s;
	}


	
	
		
}
