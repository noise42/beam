package BEAM2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class OldSearch {
	// Cerca args[0](benchmark) in args[1](fasta)

	private static void readPSSMfromBenchmark(String path, MotifHandler mh){
		FileReader f = null;

		try{
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader b = new BufferedReader(f);
			String s=b.readLine();

			while (!(s.isEmpty())){

				mh.addMotif(s.split("\t")[1], "", s.split("\t")[0], 0, s.split("\t")[0].length());
				s=b.readLine();
			}

		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


	public static void main(String[] args) {
		//cerca args[0] (PSSM) in args[1] (query)


		double[][]mbr=new double[83][83];
		BEARManager.readMatrix("./src/matriceRidondanza50",mbr); //riempie mbr
	
		System.out.println("Leggo sequenze query...");

		ArrayList<Motif> querySequences=new ArrayList<Motif>();
		
		long startTime = System.currentTimeMillis();
		IO.readInput(args[1],querySequences); //riempie querySequences (i fasta in cui cercare)
		//TODO 2
		System.out.println("elapsed Time(ms):\t" + (System.currentTimeMillis() - startTime) );
		IO.scream(+1);


		MotifHandler mh = new MotifHandler();
		double min = 0.0;
		if(args.length >= 3){
			min = Double.parseDouble(args[2]);
		}
		System.out.println("Leggo file di benchmark...");
		readPSSMfromBenchmark(args[0], mh);

		System.out.println("imposto lunghezza motivo...");
		mh.setMotifWidth();		

		System.out.println("Carico PSSM...");
		PSSM.computePSSM(mh, new ArrayList<Motif>());
		//TODO analizzare

		System.out.println("Inizio ricerca su sequenze query");
		//
		double bestMatch = 0.0;
		int counter = 0;
		for(Motif m: querySequences){
			bestMatch=0.0;
			//bisogna ciclare su le sequenze di querySequences (e per ogni sequenza su ogni subopt)
			//ciclo su sequenze
			if (m.getSequence().length() >= mh.getMotifWidth()){
				bestMatch=MotifUtilities.computeScoreVsPSSM(mh.getPSSM(), m, mbr, mh.getMotifWidth());
				//TODO 3
				//bestMatch = BEARManager.searchMotifUsingPSSMForSearch(mh.getPSSM(), m.getSequenceList(), mbr, mh.getMotifWidth());
				if (bestMatch > min){
					System.out.println(">" + m.getName());
					System.out.println(m.extractMotifFromSequence() + "\t" + m.getMotifStart() + "\t" + m.getMotifEnd() + "\t" + bestMatch);
					counter++;
				}
			}
		}
		System.out.println("\n#Ricerca completata");
		System.out.println("#numero seq query: " + counter);
		System.out.println("#su: " + querySequences.size() + " sequenze.");


		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.print("elapsed Time(ms):\t" + elapsedTime);
	}
}
