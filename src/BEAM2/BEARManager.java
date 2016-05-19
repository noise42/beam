package BEAM2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
//import java.lang.reflect.Array;
import java.util.ArrayList;

//in questa classe ci devono essere le funzioni che gestiscono la matrice e l'allineamento
public class BEARManager {

	//	static double substitutionScore(char a, char b,double[][] m){ // j k presenti!
	//		//String name="abcdefghi=jklmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";
	//		String name="abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";
	//		return m[name.indexOf(a)][name.indexOf(b)];
	//	}

	static double substitutionScore(char a, char b,double[][] m){
		//MBR
		String name="abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";
		return m[name.indexOf(a)][name.indexOf(b)];
	}

	static double substitutionScoreWithMask(char a, char b,double[][] m){
		if (a == '1' || b == '1'){
			return Double.NEGATIVE_INFINITY;
		}else{
			return substitutionScore(a, b, m);
		}
	}

	static double substitutionScore3(char a, char b,double[][] m){
		//Grezzo
		ArrayList<String> name = new ArrayList<String>();
		//Stem
		name.add("abcde");
		name.add("fghi");
		name.add("=");
		//Loop
		name.add("jklmnopqr");
		name.add("stuvwxyz");
		name.add("^");
		//Int Loop
		name.add("!\"#$%23456");
		name.add("&\'()7890");
		name.add("+>");

		//Bulge
		name.add("[]");
		//Stem B
		name.add("ABCDE");
		name.add("FGHI");
		name.add("J");
		//Int Loop B
		name.add("KLMNYZ~?");
		name.add("OPQRS_|/\\");
		name.add("TUVWYZ@");
		//Bulge B
		name.add("{}");
		//Branch
		name.add(":");
		if (a == b){
			return 1;
		}else{
			for (String classe: name){
				if (classe.indexOf(a)>=0 && classe.indexOf(b)>=0){
					return 0.5;
				}
			}
		}
		return -.1;
	}

	static double substitutionScoreInternal(char a, char b,double[][] m){
		//punti solo per gli internal loop
		String name = "!\"#$%23456&\'()7890+>[]{}KLMNYZ~?OPQRS_|/\\TUVWYZ@";

		if (name.indexOf(a)>=0 && name.indexOf(b)>=0){
			return 1;
		}else{
			return 0;
		}	
	}



	static double[][] computeScoringMatrix(double[][] mat){
		double[][] scoringMatrix = new double[83][83];
		for(int i=0; i< mat.length; i++){
			for(int j=0; j<mat.length; j++){
				scoringMatrix[i][j] = mat[i][j] - .5*(mat[i][i] + mat[j][j]);
			}
		}
		return scoringMatrix;
	}

	//	static int searchInitialMotif_(String motif,String sequence,double[][] m,int motifLength){
	//
	//		double score=0.0;
	//		int maxIndex=0;
	//
	//		for(int i=0; i< sequence.length() - motifLength; i++){
	//			double tmp=0.0;
	//			for(int j=0;j<motifLength;j++){
	//				//System.out.println(motif.charAt(j)+"\t"+sequence.charAt(j+i));
	//				tmp+=substitutionScore(motif.charAt(j),sequence.charAt(j+i), m);
	//			}
	//			if(tmp>score){
	//				score=tmp;
	//				maxIndex=i;
	//			}
	//
	//		}
	//		//System.out.println(maxIndex + "\t" + motifLength + "\t" + sequence.length());
	//
	//		return maxIndex;
	//	}

	static int searchInitialMotif(String motif,Motif m,double[][] mbr,int motifLength){

		int start=m.getMotifStart(); //backup

		double score=-99999.0;
		int maxIndex=0;
		String sequence = m.getSequence();

		for(int i=0; i< sequence.length() - motifLength; i++){
			double tmp=0.0;
			for(int j=0;j<motifLength;j++){
				//System.out.println(motif.charAt(j)+"\t"+sequence.charAt(j+i));
				tmp+=substitutionScoreWithMask(motif.charAt(j),sequence.charAt(j+i), mbr);
			}
			if(tmp>score){
				score=tmp;
				maxIndex=i;
			}

		}
		//System.out.println(maxIndex + "\t" + motifLength + "\t" + sequence.length());
		//Controllo obsoleto? Va sostituito con un controllo di fattibilita': c'e' almeno una regione lunga motifWidth non coperta?
		//TODO non vede bene quando il motifEnd e masked
		for(int i=maxIndex; i<=(maxIndex+motifLength); i++){
			if(m.isMasked(i)){
				return start;
			}
		}
		//System.out.println( m.printMask() + "\tstart: " + maxIndex + "+" + motifLength);
		return maxIndex;
	}

	static int searchInitialMotif(String motif,Motif m,double[][] mbr,int motifLength, boolean firstTime){


		double score=-99999.0;
		int maxIndex=0;
		String sequence = m.getSequence();
		int startFetch = 0;
		boolean search=true;
		while(search){

			//parti da startFetch, esci con uno score e un maxIndex
			for(int i=startFetch; i< sequence.length() - motifLength; i++){
				double tmp=0.0;
				for(int j=0;j<motifLength;j++){
					//System.out.println(motif.charAt(j)+"\t"+sequence.charAt(j+i));
					tmp+=substitutionScoreWithMask(motif.charAt(j),sequence.charAt(j+i), mbr);
				}
				if(tmp>score){
					score=tmp;
					maxIndex=i;
				}
			}
			//ora controlla se il maxIndex ha senso (non � mascherato)
			search = false;
			for(int i=maxIndex; i<=(maxIndex+motifLength); i++){
				if(m.isMasked(i)){
					startFetch++;
					//ACCROCCATISSIMO, sposta in maniera tale da ripetere la ricerca forzando lo spostamento
					search=true;
					break;
				}
			}
			//se sei uscito con un break search � di nuovo true e il punto di partenza � cambiato.
			if (sequence.length() - startFetch < motifLength){
				return -1;
			} //SE non esistono sequenze non coperte ritorna -1, gestire.
		}
		//Se esci dal while intatto con niente mascherato hai chiuso.
		//System.out.println( m.printMask() + "\tstart: " + maxIndex + "+" + motifLength);
		return maxIndex;
	}

	static int searchMotifUsingPSSM(PSSM pssm,Motif m ,double[][] mbr, int motifLength){
		int start = m.getMotifStart();
		double score=-9999.0;
		int maxIndex=0;
		String sequence_ = m.getSequence();
		if (sequence_.length()-motifLength < 0){
			return -1;
		}
		String sequence = convertUsingMasks(sequence_, m.getMask());

		for(int i=0;i<sequence.length()-motifLength;i++){
			//per ogni carattere delle sequenza, usalo come start
			double tmp=0.0;
			for(int j=0;j < motifLength;j++){
				double tmp2=0.0;
				double max2=-9999.0;
				for(PSSMCell pos:pssm.getMatrix().get(j)){

					tmp2=substitutionScoreWithMask(pos.name, sequence.charAt(j+i), mbr)*pos.occurrence;
					if(tmp2>max2){
						max2=tmp2; //prendo il massimo di ogni colonna
					}
				}
				tmp+=max2;

			}
			if(tmp>score){
				score=tmp;
				maxIndex=i;
			}
		}
		//		System.out.println(maxIndex + "\t" + motifLength + "\t" + sequence.length());
		for(int i=maxIndex; i<=(maxIndex+motifLength); i++){
			if(m.isMasked(i)){

				return start;
			}
		}
		if (score < -10000){
			System.err.println(score);
			return -1;
		}

		return maxIndex;
	}

	static int searchMotifUsingPSSM(PSSM pssm,Motif m ,double[][] mbr, int motifLength, boolean firstTime){
		double score=-9999.0;
		int maxIndex=0;
		String sequence_ = m.getSequence();
		if (sequence_.length()-motifLength < 0){
			return -1;
		}
		boolean search=true;
		int startFetch=0;

		String sequence = convertUsingMasks(sequence_, m.getMask());

		while(search){
			for(int i=startFetch;i<sequence.length()-motifLength;i++){
				//per ogni carattere delle sequenza, usalo come start
				double tmp=0.0;
				for(int j=0;j < motifLength;j++){
					double tmp2=0.0;
					double max2=-9999.0;
					for(PSSMCell pos:pssm.getMatrix().get(j)){

						tmp2=substitutionScoreWithMask(pos.name, sequence.charAt(j+i), mbr)*pos.occurrence;
						if(tmp2>max2){
							max2=tmp2; //prendo il massimo di ogni colonna
						}
					}
					tmp+=max2;

				}
				if(tmp>score){
					score=tmp;
					maxIndex=i;
				}
			}
			search=false;
			//come in searchInitialMotif
			for(int i=maxIndex; i<=(maxIndex+motifLength); i++){
				if(m.isMasked(i)){
					startFetch++;
					search=true;
					break;
				}
			}
			if (sequence.length()-startFetch < motifLength){
				return -1;
			}
		}
		return maxIndex;
	}

	static public String convertUsingMasks(String seq, ArrayList<int[]> mask){
		//converte una sequenza di un oggetto Motif in una sequenza con le zone mascherate in "1"
		char[] chars = seq.toCharArray();
		//TODO verificare se funziona
		for(int[] region: mask ){
			for (int i=region[0]; i<region[1]; i++){
				chars[i]='1';
				//stamparsi le regioni e in fondo la sequenza con gli 1 per verificare!
			}
			//System.out.println(region[0] + "\t" + region[1] + "\t\t");

		}
		//System.out.println(String.valueOf(chars));

		return String.valueOf(chars);
	}


	static String searchMotifUsingPSSMForSearch(PSSM pssm,ArrayList<String> sequenceList ,double[][] m, int motifLength){

		double score=-9999.0;
		double maxAll=0.0;
		String best = "";
		int subCount = 0;
		for(String sequence: sequenceList){
			//			if (sequence.length()-motifLength < 0){
			//				return -1;
			//			}

			for(int i=0;i<sequence.length()-motifLength;i++){
				//per ogni carattere delle sequenza, usalo come start
				double tmp=0.0;
				for(int j=0;j < motifLength;j++){
					double tmp2=0.0;
					double max2=-9999.0;
					for(PSSMCell pos:pssm.getMatrix().get(j)){

						tmp2=substitutionScore(pos.name, sequence.charAt(j+i), m)*pos.occurrence;
						if(tmp2>max2){
							max2=tmp2; //prendo il massimo di ogni colonna
						}
					}
					tmp+=max2;

				}
				if(tmp>score){
					//se la finestra cercata ora ha un punteggio piu' alto delle precedenti segna
					score=tmp;
					if(tmp>maxAll){
						maxAll=score;
						best = sequence.substring(i, i+motifLength) + "\t" + score + "\tsubopt no.:" + subCount ;
					}
				}
			}
			subCount++;
		}
		//		System.out.println(maxIndex + "\t" + motifLength + "\t" + sequence.length());


		return best;
	}

	static String searchMotifUsingPSSMForSearchWithMin(PSSM pssm,ArrayList<String> sequenceList,double[][] m, int motifLength, double min){
		//TODO vedere se serve implementare le maschere!!!
		double score=-9999.0;
		double maxAll=0.0;
		String best = "";
		int subCount = 0;
		for(String sequence: sequenceList){
			//			if (sequence.length()-motifLength < 0){
			//				return -1;
			//			}

			for(int i=0;i<sequence.length()-motifLength;i++){
				//per ogni carattere delle sequenza, usalo come start
				double tmp=0.0;
				for(int j=0;j < motifLength;j++){
					double tmp2=0.0;
					double max2=-9999.0;
					for(PSSMCell pos:pssm.getMatrix().get(j)){

						tmp2=substitutionScore(pos.name, sequence.charAt(j+i), m)*pos.occurrence;
						if(tmp2>max2){
							max2=tmp2; //prendo il massimo di ogni colonna
						}
					}
					tmp+=max2;

				}
				if(tmp>score){
					//se la finestra cercata ora ha un punteggio piu' alto delle precedenti segna
					score=tmp;
					if(tmp>maxAll){
						maxAll=score;
						best = sequence.substring(i, i+motifLength) + "\t" + score + "\tsubopt no.:" + subCount ;
					}
				}
			}
			subCount++;
		}
		//		System.out.println(maxIndex + "\t" + motifLength + "\t" + sequence.length());


		return best;
	}



	static void readMatrix(String path, double[][] m){
		FileReader f = null;

		try{
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader b = new BufferedReader(f);
			String s;
			int i=0;

			while ((s=b.readLine())!=null){
				String[] field=s.split("\t");
				for(int j = 0 ; j < field.length; j++){


					m[i][j]=Double.parseDouble(field[j]);
				}
				i++;
			}

		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	static void readVector(String path, double[] v){
		FileReader f = null;

		try{
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader b = new BufferedReader(f);
			String s;
			int i=0;

			while ((s=b.readLine())!=null){

				v[i]=Double.parseDouble(s);

				i++;
			}

		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
