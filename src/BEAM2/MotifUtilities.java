package BEAM2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;
import org.apache.commons.math3.distribution.*;

public class MotifUtilities {

	public static void computePartials(PSSM pssm, MotifHandler mh, double[][] matrix, int motifWidth){
		//input PSSM,mh. Il parziale è lo score della sequenza con la pssm, passi successivi: tenere solo quelli vicino alla media
		for (Motif m: mh.getListMotif()){
			//per ogni sequence in mh
			double tmp=0.0;
			for(int j=0; j < motifWidth;j++){
				//				scorri ogni colonna e calcola lo score come la media pesata per le occorrenze 

				for(PSSMCell pos:pssm.getMatrix().get(j)){
					tmp+=BEARManager.substitutionScore(pos.name, m.extractMotifFromSequence().charAt(j), matrix)*pos.occurrence;
				}
			}
			double partial_ = tmp;
			m.setPartial(partial_);
		}
	}

	public static double computeScoreVsPSSM(PSSM pssm, Motif m, double[][] matrix, int motifWidth){
		//input PSSM,mh. Il parziale è lo score della sequenza con la pssm, passi successivi: tenere solo quelli vicino alla media
		int idx=0;
		double max=-9999.0;
				
		for(int i=0; i < (m.getSequence().length()-motifWidth+1); i++ ){
			double tmp=0.0;
			for(int j=0; j < motifWidth;j++){
				//				scorri ogni colonna e calcola lo score come la media pesata per le occorrenze 

				for(PSSMCell pos:pssm.getMatrix().get(j)){
					tmp+=BEARManager.substitutionScore(pos.name, m.getSequence().charAt(j+i), matrix)*pos.occurrence;
					//TODO poi includere bonus sequenza
				}
			}
			//tmp2 contiene il punteggio della finestra
			if(tmp > max){
				max=tmp;
				idx=i;
			}
		}
		m.setMotifStart(idx, motifWidth);
		m.setMotifEnd(idx+motifWidth, motifWidth);
		double partial_ = max;
		return partial_;

	}

	static public double max(double[] array){
		double max = -999.0;
		for(int i=0;i<array.length; i++){
			if (array[i] > max){
				max = array[i];
			}
		}
		return max;
	}

	public static void clean(MotifHandler mh, double thr){
		double mediaPart = 0.0;
		for(Motif m: mh.getListMotif()){ //questa parte puo' essere integrata al computePartials e passato direttamente la media qua
			mediaPart += m.getPartial();
		}
		mediaPart /= mh.cardinality();
		//Ho ottenuto la media dei parziali, ora identifico quelli con punteggio meno del thr% della media
		System.out.println("mean Par:" + mediaPart);
		for(int i=(mh.cardinality()-1); i >=0 ; i--){
			if (mh.getObjectMotif(i).getPartial() < thr*mediaPart){
				//System.err.println(mh.getObjectMotif(i).getPartial());
				mh.removeMotif(i);
			}else{
				//System.out.println(mh.getObjectMotif(i).getPartial());
			}
		}
	}

	static public int maxIndex(double[] array){
		double max = -999.0;
		int maxI = 0;
		for(int i=0;i<array.length; i++){
			if (array[i] >= max){
				max = array[i];
				maxI = i;
			}
		}
		return maxI;
	}

	public static BigDecimal truncateDecimal(double x,int numberofDecimals)
	{
		if ( x > 0) {
			return new BigDecimal(String.valueOf(x)).setScale(numberofDecimals, BigDecimal.ROUND_FLOOR);
		} else {
			return new BigDecimal(String.valueOf(x)).setScale(numberofDecimals, BigDecimal.ROUND_CEILING);
		}
	}

	public static void setBranchingZero(double[][] mbr){
		mbr[48][48]=0.0;
	}

	public static void reshapeInput(ArrayList<Motif> inputSequences) {
		//controlla inputSequences se ci sono duplicati causati da boh. In caso toglili
		//Da lanciare una volta a inizio run
		ArrayList<String> list = new ArrayList<String>();
		ArrayList<Motif> toBeRem = new ArrayList<Motif>();
		for (Motif m: inputSequences){
			if (list.contains(m.getName()) ){
				toBeRem.add(m);
			}else{
				list.add(m.getName());
			}
		}
		for (Motif m: toBeRem){
			inputSequences.remove(m);
		}

	}

	public static double computeMeanLength(ArrayList<Motif> inputSequences) {
		// fornisce lunghezza media delle sequenze in inputSequences, per determinare LSBin
		double mean=0.0;
		for(Motif m: inputSequences){
			//TODO in caso di subopt gestire questa cosa
			mean+=m.getSequence().length();
		}
		return mean/inputSequences.size();
	}

	public static double computeMeanStructureContent(
			ArrayList<Motif> inputSequences) {
		// fornisce struttura media delle sequenze in inputSequences, per determinare LSBin
		double mean=0.0;
		for(Motif m:inputSequences){
			//TODO in caso di subopt gestire questa cosa
			mean+=assessStructureContent(m.getSequence());
		}
		return mean/inputSequences.size();
	}

	private static double assessStructureContent(String bear){
		double sc=0.0;
		String stemClass = "abcdefghi=ABCDEFGHIJ";
		for (int i=0; i<bear.length();i++){
			if(stemClass.contains(Character.toString(bear.charAt(i)) )  ){
				sc+=1.0;
			}
		}
		return sc/bear.length();
	}

	public static double computePvalue(String gaussFile, double meanPart) {
		// dato file Gauss + media parziali calcola zscore e poi pvalue
		double pvalue=-2;
		double mean=0.0;
		double std=0.0;
		
		FileReader f = null;
		String s="";
		
		try{
			//apro file
			f = new FileReader(gaussFile);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader reader = new BufferedReader(f);
			while ((s=reader.readLine())!=null){
				if (s.startsWith("mean")){
					mean=Double.parseDouble(s.split("\t")[1]);
				}else if(s.startsWith("var")){
					std = Math.sqrt(Double.parseDouble(s.split("\t")[1]));
				}
			}
		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		NormalDistribution gauss= new NormalDistribution(mean, std);
		pvalue= 1-gauss.cumulativeProbability(meanPart);

		


		
		return pvalue;
	}


}
