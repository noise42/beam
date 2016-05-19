package BEAM2;

import java.util.ArrayList;
//import java.util.Iterator;

public class MotifHandler{
	private ArrayList<Motif> motifList;
	private PSSM pssm, qBEARPssm;
	private int motifWidth;
	private int motifWidthPrev;
	private int shortestSequenceLength;
	private int shortestSequenceLengthPrev;
	//lunghezza sequenza piu' corta.

	public MotifHandler(){
		this.motifList=new ArrayList<Motif>();
		this.motifWidth=0;
		this.motifWidthPrev=0;
		this.pssm=new PSSM();
		this.qBEARPssm = new PSSM();
	}

	public void MotifHandlerClone(MotifHandler copy){
		this.motifList = new ArrayList<Motif>();
		for (Motif m : copy.motifList){

			Motif addM = new Motif(m);
			//			System.out.println(m.getName() +"\t" +addM.getIndex() + "\t precedente: " +m.getIndex());

			this.motifList.add(addM);
		}
		this.motifWidth=copy.motifWidth;
		this.motifWidthPrev=copy.motifWidthPrev;
		this.pssm=copy.pssm;
		this.qBEARPssm = copy.qBEARPssm;
		this.shortestSequenceLength=copy.shortestSequenceLength;
		this.shortestSequenceLengthPrev=copy.shortestSequenceLengthPrev;
	}

	public void addMotif(String name, String nuc, String seq, int s, int e){
		Motif tmp=new Motif(name, nuc, seq,s,e);
		this.motifList.add(tmp);
	}

	public void addMotif(Motif m){
		this.motifList.add(m);
	}

	public void addMotif(Motif m, int index){
		this.motifList.add(index,m);
	}
	public void removeMotif(int index){
		this.motifList.remove(index);
	}	

	public void setMotifWidth() {
		this.motifWidthPrev=this.motifWidth;
		this.motifWidth =this.motifList.get(0).getMotifEnd() - this.motifList.get(0).getMotifStart();
	}

	public void setMotifWidth(int motifWidth) {
		this.motifWidthPrev=this.motifWidth;
		this.motifWidth = motifWidth;
	}

	public int getMotifWidth() {
		return this.motifWidth;
	}

	public int getMotifWidthPrev() {
		return this.motifWidthPrev;
	}

	public Motif getObjectMotif(int index){
		return this.motifList.get(index);
	}

	public ArrayList<Motif> getListMotif(){
		return this.motifList;
	}

	public String getSequenceMotif(int i){
		return this.motifList.get(i).extractMotifFromSequence();
	}

	public int cardinality(){
		return this.motifList.size();
	}

	public String printPSSM(){
		return this.pssm.toString(this.motifWidth);
	}


	public String printQBEARPSSM(){
		return this.qBEARPssm.toString(this.motifWidth);
	}

	public double getScore(){
		return this.pssm.getScore();
	}

	static public void computeScore(MotifHandler mh){

		mh.getPSSM().setScore(PSSM.ScorePartials(mh));

	}

	public PSSM getPSSM() {
		return pssm;
	}

	public PSSM getQBEARPSSM() {
		return qBEARPssm;
	}

	@Override
	public String toString(){
		String motifString = "";
		if (Debug.ON){
			for (int i = 0; i < cardinality(); i++) {
				motifString += getSequenceMotif(i)+"\t" + getObjectMotif(i).getName() + 
						"$"+ (getObjectMotif(i).getIndex()+1) +  "\t" +  getObjectMotif(i).getSequence().length() + 
						"\t" + getObjectMotif(i).getMotifStart()+ "\t" + getObjectMotif(i).getMotifEnd() + 
						"\t" +  MotifUtilities.truncateDecimal(getObjectMotif(i).getPartial(), 4) + 
						"\t" + getObjectMotif(i).printMask() +
						"\n";
			}
		}else{
			for (int i = 0; i < cardinality(); i++) {
				motifString += getSequenceMotif(i)+"\t" + getObjectMotif(i).getName() + 
						"$"+ (getObjectMotif(i).getIndex()+1) +  "\t" +  getObjectMotif(i).getSequence().length() + 
						"\t" + getObjectMotif(i).getMotifStart()+ "\t" + getObjectMotif(i).getMotifEnd() + 
						"\t" +  MotifUtilities.truncateDecimal(getObjectMotif(i).getPartial(), 4) + 
						"\n";
			}
		}
		return motifString;
	}

	public String toString(boolean webLogo, boolean qBEARFlag){
		//in questa versione con due parametri, (webLogo) stampa il fasta, qBEARFlag converte in qBEAR
		String motifString = "";
		if(qBEARFlag){
			for (int i = 0; i < cardinality(); i++) {
				motifString += ">" + getObjectMotif(i).getName() +"\n";
				for (int j = 0; j < getSequenceMotif(i).length() ; j++){
					motifString += qBEAR.identifyClass(getSequenceMotif(i).charAt(j));
				}
				motifString += "\n";
			}
		}else{
			for (int i = 0; i < cardinality(); i++) {
				motifString += ">" + getObjectMotif(i).getName() +"\n" + getSequenceMotif(i)+"\n";
			}
		}
		return motifString;
	}

	public String toStringpBear(boolean webLogo){
		//in questa versione con due parametri, (webLogo) stampa il fasta pBear
		String motifString = "";

		for (int i = 0; i < cardinality(); i++) {
			motifString += ">" + getObjectMotif(i).getName() +"\n";
			for (int j = 0; j < getSequenceMotif(i).length() ; j++){
				//TODO motifString += qBEAR.PartialIdentifyClass(getSequenceMotif(i).charAt(j));
			}
			motifString += "\n";
		}

		for (int i = 0; i < cardinality(); i++) {
			motifString += ">" + getObjectMotif(i).getName() +"\n" + getSequenceMotif(i)+"\n";
		}

		return motifString;
	}

	public String toStringTestRegex(boolean webLogo, boolean qBEARFlag){
		//in questa versione con due parametri, (webLogo) stampa il fasta, qBEARFlag converte in qBEAR
		String motifString = "";
		if(qBEARFlag){
			for (int i = 0; i < cardinality(); i++) {
				motifString += getObjectMotif(i).getName() +"\t";
				motifString += getObjectMotif(i).getMotifStart() + "\t" + getObjectMotif(i).getMotifEnd();

				motifString += "\n";
			}
		}else{
			for (int i = 0; i < cardinality(); i++) {
				motifString += ">" + getObjectMotif(i).getName() +"\n" + getSequenceMotif(i)+"\n";
			}
		}
		return motifString;
	}

	public void setShortestSequenceLength(int shortestLength) {
		this.shortestSequenceLength = shortestLength;
	}

	public int getShortestSequenceLength() {
		return this.shortestSequenceLength;
	}

	public void setShortestSequenceLengthPrev(int shortestLength) {
		this.shortestSequenceLengthPrev = shortestLength;
	}

	public int getShortestSequenceLengthPrev() {
		return this.shortestSequenceLengthPrev;
	}

	public String toString2() {
		String motifString = "";
		for (int i = 0; i < cardinality(); i++) {
			motifString += getSequenceMotif(i);
			if (i != cardinality()-1){
				motifString += "\n";
			}
		}
		return motifString;
	}

	public String otherMatches() {
		// TODO Cerca in tutte le sequenze match del motivo e le riporta
		return "";
	}

	public String toStringSeq() {
		//stampa il motifhandler di sequenza
		String motifString = "";
		int i=0;
		try{
			for (i = 0; i < cardinality(); i++) {
				if(getObjectMotif(i).getNucleotides().equals("") == false){
					motifString += getObjectMotif(i).getNucleotides().substring(getObjectMotif(i).getMotifStart(), getObjectMotif(i).getMotifEnd())+"\t" + getObjectMotif(i).getName() + 
							"$"+ (getObjectMotif(i).getIndex()+1) +  "\tsu" +  getObjectMotif(i).getSequenceList().size() + 
							"\t" + getObjectMotif(i).getMotifStart()+ "\t" + getObjectMotif(i).getMotifEnd() + 
							"\t" +  MotifUtilities.truncateDecimal(getObjectMotif(i).getPartial(), 4) + 
							"\n";
				}
			}
		}catch(StringIndexOutOfBoundsException e){
			System.err.println("\nerrore nella sequenza numero\t" + i + "\n");
			e.printStackTrace();
		}
		return motifString;
	}

	public double getMinPartial() {
		//restituisce il minimo valore tra tutti i parziali di un mh
		double min = 999999.9;
		for (Motif m: this.motifList){
			if(m.getPartial() < min){
				min= m.getPartial();
			}
		}
		return min;
	}

	public double getMaxPartial() {
		//restituisce il max valore tra tutti i parziali di un mh
		double max = 0;
		for (Motif m: this.motifList){
			if(m.getPartial() > max){
				max= m.getPartial();
			}
		}
		return max;
	}

	public void adjustEndIndexes(int motifWidth) {
		// mette apposto gli end indexes che ogni tanto paiono scomparire
		for(Motif m:this.motifList){
			m.setMotifEnd(m.getMotifStart()+ motifWidth);
		}
		
	}





}
