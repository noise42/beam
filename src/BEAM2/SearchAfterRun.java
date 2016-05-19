package BEAM2;

import java.util.ArrayList;

public class SearchAfterRun {
	// Cerca benchmark in fasta(background + input)

	


	public static int searchAftRun(MotifHandler mh, String backgroundFile, String benchSearchFile, ArrayList<Motif> inputSequences, MotifHandler searchMh) {
		//cerca mh (PSSM) in backgroundFile+input (query)
		//restituisce grandezza background (serve per l'MCC)

		long startTime = System.currentTimeMillis();

		//DONE 1 - con questa modifica la lettura di mbr passa da 100ms a 3ms -> su 10^6 calls
		// risparmia ~ 1g 
		double[][]mbr=MBR.getmbr();
		//MotifHandler searchMh = new MotifHandler();
	
		//cose per il calcolo del pvalue
		double gaussMean=0.0, gaussVar =0.0;
		int gaussCounter = 0;
		
		//Qua Query = TARGET
		ArrayList<Motif> querySequences=new ArrayList<Motif>();
		//riempie querySequences (i fasta in cui cercare)
		//Negatives
		boolean bgFlag = true; //formale, per readability. E perch� java � na merda e non supporta i parametri opzionali.
		IO.readInputForSearch(backgroundFile,querySequences, bgFlag); 
		int bgSize = querySequences.size();
		//Positives
		for(Motif m: inputSequences){
			querySequences.add(m);

		}
		//querySequences contiene bg + input rimanente
		
		//per sicurezza setto di nuovo il parametro width del motivo (dovrebbe gi� esserci)
		mh.setMotifWidth();		

		//genero la PSSM (come sopra dovrebbe gi� esserci)
		PSSM.computePSSM(mh, new ArrayList<Motif>());
		
		MotifUtilities.computePartials(searchMh.getPSSM(), searchMh, mbr, mh.getMotifWidth());
		//leggo il max partial che c'� nel bm (sar� la threshold con cui cercare)
		double min = mh.getMinPartial();
		double max = mh.getMaxPartial();
		
		//double thr = .9*max;
		double thr = min;


		//Inizio ricerca su sequenze query
		//Nel mentre mi segno i parziali che escono dai background per avere gaussian dei parziali
		//che servira' per avere Z-score della media del modello, che -> p-value
		double bestMatch = 0.0;
		int counter = 0;
		for(Motif m: querySequences){
			bestMatch=0.0;
			
			//bisogna ciclare su le sequenze di querySequences (e per ogni sequenza su ogni subopt)
			//ciclo su sequenze
			if (m.getSequence().length() >= mh.getMotifWidth()){
				bestMatch=MotifUtilities.computeScoreVsPSSM(mh.getPSSM(), m, mbr, mh.getMotifWidth());
				
				//SEZIONE per gaussian modeling
				if (m.getName().endsWith("_bg")){
					gaussMean+=bestMatch;
					gaussVar+=bestMatch*bestMatch;
					gaussCounter+=1;
				}
				
				if (bestMatch >= thr){
					//se � accettabile a questo punto la devo integrare all'mh
					m.setPartial(bestMatch);
					mh.addMotif(m); 
					//NOTARE che start ed end sono gi� inclusi nell'oggetto Motif (sono settati da computeScoreEtc...)
					searchMh.addMotif(m);
					//System.out.println(">" + m.getName());
					//System.out.println(m.extractMotifFromSequence() + "\t" + m.getMotifStart() + "\t" + m.getMotifEnd() + "\t" + bestMatch);
					counter++;
				}
			}
		}
		System.out.println("\n#Ricerca completata");
		System.out.println("#numero seq query: " + counter);
		System.out.println("#su: " + querySequences.size() + " sequenze.");

		System.out.println(searchMh.cardinality());
		
		ArrayList<Motif> dump = new ArrayList<Motif>();
		//in dump ci vanno le sequenze che causano errori, magari pu� servire
		if( counter != 0){
			PSSM.computePSSM(searchMh, dump);
			IO.writeBenchmark(benchSearchFile, searchMh, true, -1, -1, -1);
		}else{
			IO.writeFailedSearchBenchmark(benchSearchFile);
		}
		//ora devo stamparla in un file
		
		long stopTime = System.currentTimeMillis();
		///*
		long elapsedTime = stopTime - startTime;
		System.out.print("Search - elapsed Time(ms):\t" + elapsedTime + "\n");
		System.out.print("min partial:\t" + min + "\n");
		System.out.print("max partial:\t" + max + "\n");
		System.out.print("dump size:\t" + dump.size() + "\n");

		//devo gestire il calcolo del pvalue
		gaussMean /= gaussCounter;
		gaussVar /= gaussCounter;
		gaussVar -= gaussMean*gaussMean;
		if(gaussCounter != 1){
			gaussVar *= gaussCounter/(gaussCounter-1); //correzione campionaria
		}
		String gaussfile=benchSearchFile+".gauss";
		IO.writeGaussian(gaussfile,gaussMean,gaussVar);
		
		//*/
		return bgSize;
	}

	
}
