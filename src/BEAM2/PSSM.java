package BEAM2;

//import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class PSSM {
	private ArrayList<ArrayList<PSSMCell>> matrix;
	private double score;

	public PSSM(){	
		this.matrix=new ArrayList<ArrayList<PSSMCell>>();
		this.score=0.0;
	}

	public void setScore(double d){
		this.score=d;
	}

	public double getScore(){
		return this.score;
	}

	public int getSize(){
		return matrix.size();
	}

	public ArrayList<ArrayList<PSSMCell>> getMatrix(){
		return this.matrix;
	}

	static public double[][] computeCounts(MotifHandler mh){
		//produce double[][] counts da mandare allo score


		double[][] counts = new double[mh.getMotifWidth()][83];
		String name="abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";
		int total=mh.cardinality();
		//inizializza counts
		for(int j=0;j<mh.getMotifWidth();j++){
			for(int i =0; i<83;i++){
				counts[j][i] = 0;
			}
		}
		//aumenta le count per la riga j del carattere presente sulla sequenza i
		for(int j=0;j<mh.getMotifWidth();j++){
			for(int i =0; i<total;i++){
				counts[j][name.indexOf(mh.getSequenceMotif(i).charAt(j))] += 1.0;
			}
		}

		return counts;
	}


	static public double[][] computeAlpha(MotifHandler mh, double[][] condProb, double[][] counts){
		//calcola pseudoconte per ogni colonna dell'allineamento, basandosi su mbr
		//alpha_ja = A sum_b( f_jb P(a|b) )    ---- P(a|b) = q_b exp(s(a,b)) ---- f_jb = c_jb /(sum_b' c_jb')
		int aConstant=0;
		double[][]alpha = new double[mh.getMotifWidth()][83];
		//riempire gli alpha con le q_b e le s(a,b) di MBR

		for(int j = 0; j < mh.getMotifWidth(); j++){ //per ogni colonna del motivo
			for(int a = 0; a < 83; a++){

				for(int b = 0; b < 83; b++){	
					alpha[j][a] += ( counts[j][b]/mh.cardinality() )*condProb[a][b]; 
					if(counts[j][b] > 0){
						aConstant += 1;
					}
				}
				alpha[j][a] *= 5*aConstant; //5R dovrebbe gestire il numero di sequenze

				//				System.out.println("A= " + aConstant);
				aConstant = 0; //reinizializza per ricalcolarlo per ogni alpha_ja (ottimizzabile)
			}
		}
		return alpha;

	}

	static public double sumAlpha(double[] alphaColumn){
		double tmp = 0.0;
		for(double a:alphaColumn) tmp += a;
		//		System.err.print(tmp);

		return tmp;
	}
	static public double sumCountAlpha(double[] alphaColumn, double[] countColumn){
		double tmp = 0.0;
		for(int a = 0; a < 83; a++) tmp += alphaColumn[a] + countColumn[a];

		//		System.err.print( " " + tmp);
		return tmp;
	}



	public static void computePSSM(MotifHandler mh, ArrayList<Motif> inputSequences){
		//calcola pssm di un allineamento locale di finestre, o motifHandler
		// in questa versione se una delle sequenze ha gli indici sballati la toglie

		mh.getPSSM().matrix=new ArrayList<ArrayList<PSSMCell>>();
		int total=mh.cardinality();

		mh.setMotifWidth();
		
		if (Debug.VERBOSE){
		System.out.println("[calculating PSSM...]");
		System.out.println("[processing " + total + " sequences...]");
		System.out.println("[motif width: "+ mh.getMotifWidth() +"]");
		}
		
		for(int j=0;j<mh.getMotifWidth();j++){
			ArrayList<PSSMCell> aL=new ArrayList<PSSMCell>();
			for(int i =0; i<total;i++){
				//TODO la prima condizione  un controllo per quando sballa gli indici, andrˆ risolto.
				if ((mh.getObjectMotif(i).getMotifEnd() - mh.getObjectMotif(i).getMotifStart()) > 0 && mh.getSequenceMotif(i).length() == mh.getMotifWidth() ){
					if(aL.isEmpty()){
						PSSMCell tmp=new PSSMCell(mh.getSequenceMotif(i).charAt(j), 1.0/total);
						aL.add(tmp);
					}else{
						boolean found=false;
						for(PSSMCell chiave:aL){
							if(chiave.name == mh.getSequenceMotif(i).charAt(j)){
								chiave.occurrence+=1.0/total;
								found=true;
							}
						}
						if(!found){
							PSSMCell tmp=new PSSMCell(mh.getSequenceMotif(i).charAt(j),1.0/total);
							aL.add(tmp);	
						}
					}
				}else{
					//Patch, se la sequenza sta per sballare, rimettila via e apposto
					
					inputSequences.add(mh.getObjectMotif(i));
					mh.getListMotif().remove(i);
					i--;
					total--;
				}
			}
			mh.getPSSM().matrix.add(aL);
		}
		

	}

	//*************PSSM NEW SCORE*******************
	//	//Information content, file PSSM_Theory ATTUALE 10/2014
	//		static public double InfoContent(MotifHandler mh, double[] weights, double[] dataPriors){
	//			double[][] counts = computeCounts(mh);
	//			double[][] modFreq = new double[mh.getMotifWidth()][83];
	//			double sum = 0.0;
	//			final double k = 200.0;
	//			for(int j=0;j<mh.getMotifWidth();j++){
	//				sum = 0.0;
	//				for(int a =0; a<83;a++){
	//					modFreq[j][a] = counts[j][a] + k/83;
	//					sum += counts[j][a] + k/83; //File PSSM_Theory
	//				}
	//				for(int a=0; a<83;a++){
	//					modFreq[j][a] /= sum;
	//				}
	//			}
	//			//ho le f'_ij
	//			double infoContent = 0.0;
	//			double infotmp=0.0;
	//			for(int j=0;j<mh.getMotifWidth();j++){
	//				for(int a =0; a<83;a++){
	//					infotmp += modFreq[j][a] * Math.log(modFreq[j][a]/dataPriors[a])*weights[a];
	//
	//				}
	//				infoContent += infotmp;
	//			}
	//			return infoContent;
	//		}

	public static void computePSSM_(MotifHandler mh){
		//calcola pssm di un allineamento locale di finestre, o motifHandler
	
	
		mh.getPSSM().matrix=new ArrayList<ArrayList<PSSMCell>>();
		int total=mh.cardinality();
	
		mh.setMotifWidth();
		
		if (Debug.VERBOSE){
		System.out.println("[calculating PSSM...]");
		System.out.println("[processing " + total + " sequences...]");
		System.out.println("[motif width: "+ mh.getMotifWidth() +"]");
		}
		
		for(int j=0;j<mh.getMotifWidth();j++){
			ArrayList<PSSMCell> aL=new ArrayList<PSSMCell>();
			for(int i =0; i<total;i++){
				if (mh.getSequenceMotif(i).length() == mh.getMotifWidth()){
					if(aL.isEmpty()){
						PSSMCell tmp=new PSSMCell(mh.getSequenceMotif(i).charAt(j), 1.0/total);
						aL.add(tmp);
					}else{
						boolean found=false;
						for(PSSMCell chiave:aL){
							if(chiave.name == mh.getSequenceMotif(i).charAt(j)){
								chiave.occurrence+=1.0/total;
								found=true;
							}
						}
						if(!found){
							PSSMCell tmp=new PSSMCell(mh.getSequenceMotif(i).charAt(j),1.0/total);
							aL.add(tmp);	
						}
					}
				}else{
					
				}
			}
			mh.getPSSM().matrix.add(aL);
		}
	}

	//Information content, file PSSM_Theory ATTUALE 10/2014 INTEGRAZIONE SEQUENZA, BACKUP DI SOPRA
	static public double InfoContent(MotifHandler mh, double[] weights, double[] dataPriors, double bonusSeq, double nThr){
		double[][] counts = computeCounts(mh);
		double bonus = 0.0;
//		if (bonusSeq != 0.0 && mh.cardinality() > 4){
//			bonus = computeBonusSeq(mh, bonusSeq, nThr);
//		}
		double[][] modFreq = new double[mh.getMotifWidth()][83];
		double sum = 0.0;
		final double k = 83*2;
		for(int j=0;j<mh.getMotifWidth();j++){
			sum = 0.0;
			for(int a =0; a<83;a++){
				modFreq[j][a] = counts[j][a] + k/83;
				sum += counts[j][a]; //File PSSM_Theory
			}
			for(int a=0; a<83;a++){
				modFreq[j][a] /= (sum+k);
			}
		}
		//ho le f'_ij
		double infoContent = 0.0;
		double infoCol=0.0;
		for(int j=0;j<mh.getMotifWidth();j++){
			for(int a =0; a<83;a++){
				infoCol += modFreq[j][a] * Math.log(modFreq[j][a]/dataPriors[a])*weights[a];

			}
			infoContent += infoCol;
		}
		return infoContent + bonus;
	}

	private static double computeBonusSeq(MotifHandler mh, double bonusSeq, double nThr) {
		//calcola le frequenze dei NUCLEOTIDI nell'allineamento 
		double[][] seqFreq = new double[mh.getMotifWidth()][4];
		String name="ACTG";
		double bonus = 0.0;
		int total=mh.cardinality();
		//inizializza counts
		for(int j=0;j<mh.getMotifWidth();j++){
			for(int i=0; i<4;i++){
				seqFreq[j][i] = 0.0;
			}
		}
		//aumenta le count per la riga j del carattere presente sulla sequenza i
		for(int j=0;j<mh.getMotifWidth();j++){//per ogni colonna controlla se c' conservazione > di 70%
			for(int i=0; i<total;i++){
				seqFreq[j][name.indexOf(mh.getObjectMotif(i).getNucleotidesMotif().toUpperCase().replace("U", "T").charAt(j)  )] += 1.0;
			}

			boolean consCol = false;
			for(int k =0; k<4;k++){
				if (seqFreq[j][k] /total > nThr){
					consCol = true;
				}
			}
			if (consCol == true){
				bonus += bonusSeq;
				consCol = false;
			}
		}
		return bonus;
	}
	
	public static BigDecimal checkBonusSeq(MotifHandler mh) {
		//calcola le frequenze dei NUCLEOTIDI nell'allineamento 
		double[][] seqFreq = new double[mh.getMotifWidth()][15];
		String name="ACGUTRYMKSWBDHVN";
		double bonus = 0.0;
		int total=mh.cardinality();
		//inizializza counts
		for(int j=0;j<mh.getMotifWidth();j++){
			for(int i=0; i<15;i++){
				seqFreq[j][i] = 0.0;
			}
		}
		//aumenta le count per la riga j del carattere presente sulla sequenza i
		for(int j=0;j<mh.getMotifWidth();j++){//per ogni colonna controlla se c' conservazione > di 70%
			for(int i=0; i<total;i++){
				seqFreq[j][name.indexOf(mh.getObjectMotif(i).getNucleotidesMotif().toUpperCase().replace("U", "T").charAt(j)  )] += 1.0;
			}

			boolean consCol = false;
			for(int k =0; k<15;k++){
				if (seqFreq[j][k] /total > 0.7){
					consCol = true;
				}
			}
			if (consCol == true){
				bonus += 1.0;
				consCol = false;
			}
		}
		return MotifUtilities.truncateDecimal(bonus/mh.getMotifWidth(), 2);
	}
	//************FINE PSSM NEW SCORE****************


	//*************SCORE CON PARZIALI*************
	public static double ScorePartials(MotifHandler mh){ //base, senza bonus seq, non utilizzato
		double score = 0.0;
		for (Motif m:mh.getListMotif()){
			score += m.getPartial();
		}
		return score;
	}

	public static double ScorePartials(MotifHandler mh, double bonusSeq, double nThr){
		//deprecated
		double score = 0.0;
		for (Motif m:mh.getListMotif()){
			score += m.getPartial();
		}
		
		if (bonusSeq!=0.0 && mh.cardinality() >= 10){
			Debug.check(1);
			score += computeBonusSeq(mh,bonusSeq, nThr);
		}
		return score;
	}
	//**********FINE SCORE CON PARZIALI********


	public class CustomComparator implements Comparator<PSSMCell> {

		@Override
		public int compare(PSSMCell arg0, PSSMCell arg1) {
			//Implementato al contrario per avere un reverse sort ;) 
			if(arg0.occurrence - arg1.occurrence > 0){
				return -1;
			}else if(arg0.occurrence - arg1.occurrence < 0){
				return 1;
			}else{
				return 0;
			}
		}
	}

	public String toString(int motifWidth){
		String out="";
		for(int j=0;j<motifWidth;j++){
			Collections.sort(this.matrix.get(j), new CustomComparator());
			for(PSSMCell chiave:this.matrix.get(j) ){
				out+=chiave.name+" : "+MotifUtilities.truncateDecimal(chiave.occurrence,2)+"\t";
			}
			out+="\n";
		}
		return out;
	}



	public static void computeQBEARPSSM(MotifHandler mh){
		mh.getQBEARPSSM().matrix=new ArrayList<ArrayList<PSSMCell>>();
		int total=mh.cardinality();
		for(int j=0;j<mh.getMotifWidth();j++){	//su ogni carattere del motivo
			ArrayList<PSSMCell> aL=new ArrayList<PSSMCell>();
			for(int i =0; i<total;i++){ //per ogni sequenza
				if(aL.isEmpty()){
					PSSMCell tmp=new PSSMCell(qBEAR.identifyClass(mh.getSequenceMotif(i).charAt(j) ), 1.0/total);
					aL.add(tmp);
				}else{
					boolean found=false;
					for(PSSMCell chiave:aL){
						//System.out.println(this.motifList.get(i).getMotifStart()+"\t"+this.motifList.get(i).getMotifEnd()+"\t"+this.motifList.get(i).getSequence());
						if(chiave.name == qBEAR.identifyClass(mh.getSequenceMotif(i).charAt(j))){
							chiave.occurrence+=1.0/total;
							found=true;
						}
					}
					if(!found){
						PSSMCell tmp=new PSSMCell(qBEAR.identifyClass(mh.getSequenceMotif(i).charAt(j)),1.0/total);
						aL.add(tmp);	
					}
				}
			}
			mh.getQBEARPSSM().matrix.add(aL);
		}
	}

}
