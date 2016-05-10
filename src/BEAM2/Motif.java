package BEAM2;

import java.util.ArrayList;

public class Motif{
	private String name;
	private String nucleotides;
	private String dotBracket;
	private String sequence;
	private ArrayList<String> sequence_;
	private int sequenceIndex;
	private int sequenceIndexPrev=0; //tiene traccia della precedente subopt (serve per l'undo)
	private int motifStart;
	private int motifStartPrev=0; //tiene traccia del precedente start (serve per l'undo)
	private int motifEnd;
	private int motifEndPrev=0; //tiene traccia del precedente end (serve per l'undo)
	private double partial;
	private double partialPrev=0;
	private ArrayList<ArrayList<int[]>> mask_; //per ogni subopt ci sar� assegnata un maschera, ovvero una coppia di posizioni da evitare
	
	
	public Motif(){
		this.sequence_ = new ArrayList<String>();
		this.sequence_.add("");
		this.mask_ = new ArrayList<ArrayList<int[]>>();
		this.setName("");
		this.setNucleotides("");
		this.setDotBracket("");
		this.setSequence("");
		this.setMotifStart(0);
		this.setMotifEnd(0);
		this.setIndex(0);
		this.setPartial(0.0);
	}

	public Motif(Motif copy){
		this.setName(copy.name);
		this.setNucleotides(copy.nucleotides);
		this.setDotBracket(copy.dotBracket);
		this.mask_ = copy.getMaskList(); //quando si copia un motivo serve anche la maschera

		this.setSequence(copy.sequence);
		this.sequence_ = new ArrayList<String>(copy.getSequenceList());
		this.setIndex(copy.sequenceIndex);
		this.setMotifStart(copy.motifStart);
		this.setMotifEnd(copy.motifEnd);
		this.setPartial(copy.partial);
	}

	public Motif(String n,String nuc, String seq,int s,int e){
		this.mask_ = new ArrayList<ArrayList<int[]>>();

		this.setName(n);
		this.setNucleotides(nuc);
		this.setDotBracket("");
		this.sequence_ = new ArrayList<String>();
		this.sequence_.add(seq);
		this.setMotifStart(s);
		this.setMotifEnd(e);
		this.setIndex(0);
		this.setPartial(0.0);

	}
	
	public Motif(String n,String nuc, String dotB, String seq,int s,int e){
		this.mask_ = new ArrayList<ArrayList<int[]>>();

		this.setName(n);
		this.setNucleotides(nuc);
		this.setDotBracket(dotB);
		this.sequence_ = new ArrayList<String>();
		this.sequence_.add(seq);
		this.setMotifStart(s);
		this.setMotifEnd(e);
		this.setIndex(0);
		this.setPartial(0.0);

	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}
	
	
	public String getNucleotides() {
		return nucleotides;
	}
	
	public String getNucleotidesMotif(){
		return nucleotides.substring(motifStart, motifEnd);
	}


	public void setNucleotides(String nucleotides) {
		this.nucleotides = nucleotides;
	}
	
	public String getDotBracket() {
		return dotBracket;
	}

	public void setDotBracket(String dotBracket) {
		this.dotBracket = dotBracket;
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
		
		//this.mask_.add(new ArrayList<int[]>());

	}

	public void setSequence(String sequence, int index) {
		//setta la sequenza all'indice index dell'arrayList
		this.sequence_.set(index, sequence);
		
		this.mask_.add(new ArrayList<int[]>());
	}



	public String getSequence() {
		//ritorna l'attuale subopt scelta
		return sequence_.get(sequenceIndex);
	}

	public String getSequence(int subopt) {
		//versione a piu' sequenze per gestire piu strutture suboptimali
		return sequence_.get(subopt);
	}

	public ArrayList<String> getSequenceList(){
		return sequence_;
	}
	
	public void setIndex(int index){
		//setta l'index della subopt scelta
		this.sequenceIndexPrev = this.sequenceIndex;
		this.sequenceIndex = index;
		//System.out.println("in caso di ripristino: " + this.sequenceIndexPrev + "\t attuale: "+ this.sequenceIndex);
	}
	public int getIndex(){
		return this.sequenceIndex;
	}
	
	public int getIndexPrev(){
		return this.sequenceIndexPrev;
	}

	public void clearStartEndIndexes(){
		this.motifEnd=0;
		this.motifEndPrev=0;
		this.motifStart=0;
		this.motifStartPrev=0;
	}

	private int getMotifWidth(){
		return this.motifEnd - this.motifStart;
	}
	
	public boolean setMotifStart(int motifStart) {
		///*
		if(isMasked(motifStart) || isMasked(motifStart + this.getMotifWidth())){ //ancora incerto per il resize (MA dovrebbe essere gestito a monte nel metodo)
			return false; //annuncio che ho provato a settare uno start non fattibile
		}//*/
		if(motifStart >= 0){
			this.motifStartPrev=this.motifStart;
			this.motifStart = motifStart;
		} else {
			System.err.println("Index Start minore di 0 !");
			System.exit(-1);
		}
		
		return true;
		
	}
	
	public boolean setMotifStart(int motifStart, int motifWidth){
		///*
		if(isMasked(motifStart) || isMasked(motifStart + motifWidth)){ //ancora incerto per il resize (MA dovrebbe essere gestito a monte nel metodo)
			this.motifStartPrev=this.motifStart;
			return false; //annuncio che ho provato a settare uno start non fattibile
		}//*/
		if(motifStart >= 0){
			this.motifStartPrev=this.motifStart;
			this.motifStart = motifStart;
		} else {
			System.err.println("Index Start minore di 0 !");
			return false;
		}
		
		return true;
		
	}

	public int getMotifStart() {
		return motifStart;
	}


	public boolean setMotifEnd(int motifEnd) {
		///*
		if(isMasked(motifEnd) || isMasked(motifEnd - this.getMotifWidth() )  ){
			//System.err.println("masked");
			return false; //annuncio che ho provato a settare un end non fattibile
		}//*/
		if(motifEnd <= this.getSequence().length()){
			this.motifEndPrev=this.motifEnd;
			this.motifEnd = motifEnd;
		}else if(motifEnd <= 0){
			System.err.println(this.name + "\t" + motifEnd);
			System.err.println("Index End minore o uguale a 0 !");
			System.exit(-1);
		}else{
			System.err.println(this.name + "\t" + this.sequence_.get(sequenceIndex).length() + "\t" + motifEnd);
			System.err.println("Index End maggiore della lunghezza massima !");
			System.err.println("Index sub: " + this.sequenceIndex);

			System.exit(-1);
		}
		
		return true;
		
	}
	public boolean setMotifEnd(int motifEnd,int motifWidth) {
		///*
		if(isMasked(motifEnd) || isMasked(motifEnd - motifWidth )  ){
			this.motifEndPrev=this.motifEnd;
			//System.err.println("masked");
			return false; //annuncio che ho provato a settare un end non fattibile
		}//*/
		if(motifEnd <= this.getSequence().length()){
			this.motifEndPrev=this.motifEnd;
			this.motifEnd = motifEnd;
		}else if(motifEnd <= 0){
			System.err.println(this.name + "\t" + motifEnd);
			System.err.println("Index End minore o uguale a 0 !");
			System.exit(-1);
		}else{
			System.err.println(this.name + "\t" + this.sequence_.get(sequenceIndex).length() + "\t" + motifEnd);
			System.err.println("Index End maggiore della lunghezza massima !");
			System.err.println("Index sub: " + this.sequenceIndex);

			System.exit(-1);
		}
		
		return true;
		
	}
	
	
	public int getMotifEnd() {
		return motifEnd;
	}

	public void setPartial(double partial_) {
		this.partialPrev = this.partial;
		this.partial = partial_;
	}	
	
	public double getPartial() {
		return this.partial;
	}

	
	public void setMotifStartUndo() {
		this.motifStart=this.motifStartPrev;
	}

	public void setMotifEndUndo() {
		this.motifEnd = this.motifEndPrev;
	}

	public void setIndexUndo() {
		this.sequenceIndex = this.sequenceIndexPrev;
	}

	public void setPartialUndo() {
		this.partial = this.partialPrev;
	}	
	
	public String extractMotifFromSequence() {
		return this.getSequence().substring(motifStart, motifEnd);
	}

	
	//MAIN method
	public void addSuboptMask(){
		this.mask_.add(new ArrayList<int[]>());		
	}
	
	public void addMask(int start, int end){
		int[] se = new int[2];
		se[0] = start;
		se[1] = end;
		//this.mask_.get(sequenceIndex).add(se);
		this.mask_.get(0).add(se);

	}
	
	public ArrayList<ArrayList<int[]>> getMaskList(){
		return this.mask_;
	}
	
	public ArrayList<int[]> getMask(){
		return this.mask_.get(0);
	}
	
	public String printMask(){
		String s = name;
		for(ArrayList<int[]> al: mask_){
			//s += "\nsub:\n";
			for(int[] msk: al){
				s += "\t[" + msk[0] + "," + msk[1]+ "]\t";
			}
		}
		
		return s;
	}
	
	@Override
	public String toString() {
		String motifString = "";
		motifString=">"+this.name + "|" + this.motifStart + "|" + this.motifEnd + "\n" + this.getSequence() + "\n";
		return motifString;
	}

	public boolean equals(Motif m){
		// if the two objects are equal in reference, they are equal
		if (this == m) {
			return true;
		} else if (m instanceof Motif) {
			if (m.getName().equals(this.getName()) && m.getSequence().equals(this.getSequence())) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
	
	public boolean isMasked(int position){
		//if position is in any of the choosen subopt masks then return true.
		if(mask_.isEmpty()){ //non fare il controllo se la maschera non esiste
			return false;
		}
		boolean isMasked_= false;
				for(int[] se: mask_.get(0)){
					if(position >= se[0] && position < se[1]){
						isMasked_=true;
						return isMasked_;
					}
				}
			return isMasked_;
	}

	public int getMotifEndPrev() {
		// TODO Auto-generated method stub
		return motifEndPrev;
	}

	public int getMotifStartPrev() {
		// TODO Auto-generated method stub
		return motifStartPrev;
	}



}
