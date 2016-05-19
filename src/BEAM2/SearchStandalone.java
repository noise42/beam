package BEAM2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class SearchStandalone {
	// Cerca args[0](benchmark) in args[1](fasta)

	public static void readPSSMfromBenchmark(String path, MotifHandler mh){
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

		long startTime = System.currentTimeMillis();

		//DONE 1 - con questa modifica la lettura di mbr passa da 100ms a 3ms -> su 10^6 calls
		// risparmia ~ 1g 
		double[][]mbr=MBR.getmbr();
		
	
		
		System.out.println("Leggo sequenze query...");

		ArrayList<Motif> querySequences=new ArrayList<Motif>();

		IO.readInputForSearch(args[1],querySequences); //riempie querySequences (i fasta in cui cercare)
		//TODO upgradabile (risparmio di qualche ora per ora)


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
		///*
		long elapsedTime = stopTime - startTime;
		System.out.print("elapsed Time(ms):\t" + elapsedTime);
		//*/
	}


	public static void computeMCC(String baseName, int mask, int bestRunIdx) {
		// TODO dal motifHandler finale lancia search su background se fornito, 
		//altrimenti su rfam - criteri da stabilire
		
	}
}
