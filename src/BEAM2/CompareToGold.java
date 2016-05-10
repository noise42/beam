package BEAM2;
import java.util.ArrayList;

public class CompareToGold {
	//questa classe fornisce un main per le analisi di BEAM. Compara la pssm finale con la string BEAR gold standard che va fornita dall'esterno

	private static double Compare(ArrayList<ArrayList<PSSMCell>> pssmMatrix, String gold, double[][] mbr){
		//Questo metodo confronta una pssm con una sequenza determinandone la vicinanza, calcolando lo score tramite MBR di un confronto pssm-seq e di uno seq-seq e facendone il rapporto

		if (gold.length() != pssmMatrix.size()){
			System.err.println("lunghezze gold standard e pssm non corrispondono\n" + "gold= "+ gold.length() +"\n" + "pssm= " + pssmMatrix.size() );
			System.exit(-1);
		}

		double score = 0.0;
		for(int i=0; i<pssmMatrix.size(); i++){
			double tmp=0.0;
			for(int j=0; j<pssmMatrix.get(i).size(); j++){
				tmp+= BEARManager.substitutionScore(pssmMatrix.get(i).get(j).name, gold.charAt(i), mbr)*pssmMatrix.get(i).get(j).occurrence;
			}
			double rowScore=tmp;
			score += rowScore;
		}



		//max
		double max = 0.0;
		for(int i=0; i<gold.length(); i++){
			max += BEARManager.substitutionScore(gold.charAt(i),gold.charAt(i),mbr);
		}
		double cmp = score/max;

		return cmp;
	}

	private static double bestCompare(PSSM pssm, String gold, double[][] mbr){
		//se gold standard e motivo finale della run sono di lunghezze diverse confronta tutta le finestre possibili e restituisci la più alta.
		if (gold.length() == pssm.getMatrix().size()){
			return Compare(pssm.getMatrix(), gold, mbr);
		}else if(gold.length() > pssm.getMatrix().size()){
			System.out.println("[looking for the best window of size " + pssm.getMatrix().size() + "...]");
			double max=0.0;
			double tmp=0.0;
			for(int i=0; i<= gold.length()-pssm.getMatrix().size(); i++){
				tmp = Compare(pssm.getMatrix(), gold.substring(i, pssm.getMatrix().size()+i), mbr);
				if (tmp>max) max = tmp;
			}
			return max;
		}else{
			System.out.println("[looking for the best window of size " + gold.length() + "...]");

			double max=0.0;
			double tmp=0.0;
			for(int i=0; i<= pssm.getMatrix().size() - gold.length(); i++){
				//Necessario convertire List in ArrayList
				ArrayList<ArrayList<PSSMCell>> pssmSubMatrix = new ArrayList<ArrayList<PSSMCell>>(pssm.getMatrix().subList(i, gold.length()+i));
				tmp = Compare(pssmSubMatrix, gold, mbr);
				
				//System.out.println(">>>tmp " + tmp + " <<<");
				
				if (tmp>max) max = tmp;
			}
			return max;
		}
	}

	public static void main(String[] args) {
		//read the gold standard (string for debug)

		String gold = args[0];


		//read the motifHandler in input file (output of run)
		MotifHandler mh = new MotifHandler();
		SCORE.readMotifHandler(args[1], mh);

		double[][] mbr = new double[83][83];
		BEARManager.readMatrix("./src/matriceRidondanza50",mbr); //riempie mbr
		
		PSSM.computePSSM(mh, new ArrayList<Motif>());

		double cmp= bestCompare(mh.getPSSM(), gold, mbr);

		System.out.println("gold standard " + gold);


		System.out.println("il motivo è corretto per il " + MotifUtilities.truncateDecimal(cmp*100, 2) + "%");
	}

}
