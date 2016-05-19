package BEAM2;

//import java.util.ArrayList;
import java.io.BufferedReader;
//import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class SCORE {

	public static void readMotifHandler(String path, MotifHandler mh){
		FileReader f = null;

		try{
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader b = new BufferedReader(f);
			String s=b.readLine();
			int i=0;
			while (!s.equals("")){ 
//				System.out.println(s);

				mh.addMotif(">"+i, "", s.split("\t")[0], 0, s.split("\t")[0].length());	
				i++;
				if((s=b.readLine()) ==null){
					break;
				}

			}

		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		if(mh.getMotifWidth() == 0){
			//Quando l'mh viene riempito da file esterno non ha settato il motif width e i motifEnd
			mh.setMotifWidth(mh.getObjectMotif(0).getSequence().length());
		}
	}


	public static void main(String[] args) {


		double[][]mbr=new double[83][83];
		double[][]condProb = new double[83][83];
		double[]dataPriors = new double[83];
		double[]alphas = new double[83];
		BEARManager.readMatrix("./src/matrixPriori", condProb);
		BEARManager.readMatrix("./src/matriceRidondanza50",mbr); //riempie mbr
		BEARManager.readVector("./src/listaPriorPesati_uniform.txt", dataPriors);
		BEARManager.readVector("./src/alpha.txt", alphas); //mettere pesi

		MotifHandler mh = new MotifHandler();

		readMotifHandler(args[0], mh);

		mh.setMotifWidth();		

		PSSM.computePSSM(mh, new ArrayList<Motif>());

//		MotifHandler.computeScore(mh, alphas, dataPriors, 0.0, 0);//genera primo punteggio

		mh.setMotifWidth();

		PSSM.computePSSM(mh, new ArrayList<Motif>());

	//	MotifHandler.computeScore(mh, alphas, dataPriors, 0.0, 0);


		System.out.println(mh.toString());
		System.out.println("------");
		System.out.println("Motif Score : "+mh.getScore());


	}

}