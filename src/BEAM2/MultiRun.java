package BEAM2;

//
//Dear maintainer:
//	
//When I wrote this code, only I and Eugenio and God 
//knew what it was. 
//
//Now, only God knows!
//
//So if you are done trying to 'optimize' 
//this routine (and failed),
//please increment the following counter
//as a warning
//to the next guy:
//	
//total_hours_wasted_here = 0
//
//Sincerely, 
//Marco "Noise" Pietrosanto
//

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class MultiRun extends BEAM2 {

	public static void main(String[] args) throws IOException {
		MotifMemory mm = new MotifMemory();
		MotifMemory tmpMm = new MotifMemory();
		SummarizedData sumData = new SummarizedData();
		//genero inputSequences
		ArrayList<Motif> inputSequences=new ArrayList<Motif>();

		int numberOfRuns = 1;
		int numberOfMasks = 1;
		String baseName = "lastrun";
		String backgroundInput = ""; //DEFAULT rfam?
		boolean USERBG=false;
		boolean WARNING_OOBIN=false; //true se USERBG false & mean length > 500
		//CommandLine Parser 
		CommandLineParser clp = new CommandLineParser(args);
		//initialized CLP, the hash is created, check with - GetValue(key), ContainsKey(key)
		if(clp.containsKey("R")){
			numberOfRuns = Integer.parseInt(clp.getValue("R"));
		}
		if(clp.containsKey("M")){
			numberOfMasks = Integer.parseInt(clp.getValue("M"));
		}
		if(clp.containsKey("f")){
			String[] tmp =  clp.getValue("f").split("/");
			baseName = tmp[tmp.length-1].split("\\.")[0]; //
			IO.readInput(clp.getValue("f"), inputSequences);
		}
		if(clp.containsKey("v")){
			Debug.setVERBOSE(true);
		}
		if(clp.containsKey("k")){
			//se esiste -k (keep) allora tieni i risultati di tutte le run
			IO.keepRuns();
		}
		if(clp.containsKey("g")){
			//se esiste -k (keep) allora tieni i risultati di tutte le run
			backgroundInput = clp.getValue("g");
			USERBG = true;
		}



		//

		String outFolder="risultati/" + baseName;
		final File dir = new File(outFolder);
		if(dir.exists()){
			boolean parent = false;
			IO.deleteFolder(outFolder, parent);
		}
		/*
			System.out.println("My scans tell me I already own an output with this name -"+baseName+ "-, want me to overwrite it?");

			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(System.in));
			String response = bufferedReader.readLine();
			//String response = console.readLine("y/n");
			if(response.toLowerCase().equals("y")){
				System.out.println("OVERWRITING...");
				IO.deletePrevious(outFolder);
			}else{
				System.out.println("CLOSING...");
				System.exit(0);
			}
		} 
		 */
		String outBench=outFolder+"/benchmark";
		String outBest = outBench + "/motifs";
		String outLogo=outFolder+"/webLogoOut";
		String outBestLogo= outLogo + "/motifs";
		String summaryName= outFolder+"/" + baseName + "_summary.txt";
		dir.mkdirs();


		new File(outBench).mkdirs();
		new File(outBest).mkdirs();
		new File(outLogo).mkdirs();
		new File(outBestLogo).mkdirs();



		int bestRunIdx=0;

		long startTimeGlobal = System.currentTimeMillis();


		if(!USERBG){
			//se non e' presente un background inserito dall'utente, prova comunque a stimarlo
			//da Rfam, con i criteri di LSbin
			String bin=sumData.computeBin();
			if (bin.charAt(0) == '4') WARNING_OOBIN=true;
			backgroundInput = IO.defaultBgCreation(bin, inputSequences.size(), outFolder);			
		}

		//Riempio i parametri iniziali per il summary (possibili solo input e bg a questo punto)
		//ma anche i parametri per il bin
		sumData.setInputFile(baseName);
		String[] bgTmp = backgroundInput.split("/");
		sumData.setBgFile(bgTmp[bgTmp.length - 1].split("\\.")[0]);
		sumData.setMl(MotifUtilities.computeMeanLength(inputSequences));
		sumData.setMsc(MotifUtilities.computeMeanStructureContent(inputSequences));

		//NOTA: sono anche gli unici fissi (gli altri cambiano a ogni run)
		
		for(int mask=1; mask<=numberOfMasks; mask++){
			long startTime = System.currentTimeMillis();

			ArrayList<Motif> clonedInput = new ArrayList<Motif>();
			mm=null;
			mm=new MotifMemory();

			tmpMm = null;
			tmpMm =new MotifMemory();

			String outMask = outBench + "/mask" + mask;
			String outMaskLogo = outLogo + "/mask" + mask;

			new File(outMask).mkdirs();
			new File(outMaskLogo).mkdirs();

			
			for(int i=1; i<=numberOfRuns; i++){

				tmpMm=run(args, i, mask, outFolder, inputSequences);
				//se lo score � pi� alto, allora metti in cloned le inputSeq (-motifhandl)
				if(mm.tryMask(tmpMm.getScore(), tmpMm.getHandlerMemory())){
					bestRunIdx = i;
					clonedInput = new ArrayList<Motif>(inputSequences.size()); //Se la maschera e' l'attuale migliore, segnati anche lo stato di input sequences
					for(Motif motif: inputSequences){
						clonedInput.add(new Motif(motif));
					}
				}

				if(i != numberOfRuns){ //altrimenti all'ultima run della maschera refilla due volte -- all'ultimo va refillato col migliore
					refillInput(tmpMm.getHandlerMemory(), inputSequences);
				}
				//outFolder � basename/
				System.out.println("mask: " + mask +"\trun: " + i + "\tcompleted");
			}
			//in bestRunIdx e' contenuto l'indice della run da copiare, in mask c'e' invece il numero di maschera






			//applyMask(mm.getHandlerMemory());
			refillInput(mm.getHandlerMemory(), clonedInput); //Applica la maschera ai Motif del migliore mh e riempi l'inputSequences associato ad esso con i motif 
			//aggiornati

			//Dopo aver deciso le migliori...
			//SEARCH -- dopo il refill cos� il background pu� essere composto
			//facilmente anche delle input refillate
			MotifHandler searchMh = new MotifHandler();

			if(!backgroundInput.equals("")){
				String backgroundFile = backgroundInput;
				String searchBenchName = outBest + "/" + baseName + "_m" + mask + "_run" + bestRunIdx + ".search.txt";
				int bgSize = SearchAfterRun.searchAftRun(mm.getHandlerMemory(), backgroundFile, searchBenchName, clonedInput, searchMh);	
				//NOTA: in questo modo se in futuro voglio tener conto per la mask anche delle
				//sequenze aggiunte con search NON posso perch� ho gia fatto il refill.
				//per� questa scelta mi serve per montare facilmente input + bg

				SearchObj sObj = IO.readSearchOutput(searchBenchName, clonedInput.size(), bgSize );
				
				
				String gaussFile = searchBenchName + ".gauss";
				double pvalue;
				pvalue=MotifUtilities.computePvalue(gaussFile, sObj.getMeanPart());
				sObj.setPvalue(pvalue);
				//Aggiornamento sumData
				sumData.setMCC(sObj.computeMCC());
				sumData.setScore(mm.getHandlerMemory().getScore());
				sumData.setConsensus(sObj.computeConsensus());
				sumData.setqbearConsensus(sObj.computeqbearConsensus());

				sumData.setInputSize(sObj.getTP() + sObj.getFN());
				sumData.setBgSize(sObj.getFP()+sObj.getTN());
				sumData.setMask(mask);
				sumData.setCoverage(sObj.getTP()*1.0 / (sObj.getTP()+ sObj.getFN()) ) ;
				sumData.setFallout(sObj.getFP()*1.0/ (sObj.getFP()+sObj.getTN()));
				sumData.setAcc((sObj.getTP() + sObj.getTN())*1.0/(sObj.getTP()+ sObj.getFN()+sObj.getFP()+sObj.getTN()));
				
				
				long elapsedTime = System.currentTimeMillis() - startTime;
				IO.appendToSummary(summaryName, sObj, sumData, elapsedTime);
			}
			
			//la mask la applico dopo il search
			if(mask != numberOfMasks){ 
				System.out.println("applying mask");
			}
			applyMask(searchMh);
			

			inputSequences = null; //svuota inputSequences e riempila copiando gli elementi di ClonedInput
			inputSequences = new ArrayList<Motif>();
			//prova con searchmh, dopo bisogner� integrare con quelli rimasti fuori
			for(Motif motif: searchMh.getListMotif()){
				if(!motif.getName().endsWith("_bg")){
					inputSequences.add(motif);
					//System.out.println(motif.getName());
				}
			}
			for(Motif motif: clonedInput){
				if(!inputSequences.contains(motif)){
					
					inputSequences.add(new Motif(motif));
				}
				/*
				else{
					System.out.println("ho gi� " + motif.getName() + ", con maschera " + motif.printMask());
				}//QUAQUA */
			}
			//System.out.println("in partenza avevo " + startingSize + ", ora ho " + inputSequences.size());

			//COPIA cartelle
			IO.moveBestRuns(baseName, mask, bestRunIdx, outBench);
			IO.moveBestRunsLogo(baseName, mask, bestRunIdx, outLogo);
			//weblogo
			Weblogo weblogo = new Weblogo();
			try {
				weblogo.generateLogo(baseName, mask, bestRunIdx, outLogo);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			SearchStandalone.computeMCC(baseName, mask, bestRunIdx);
			IO.appendToBenchmark(baseName,mask,bestRunIdx);

			//Cancella cartelle inutili se non � stato richiesto di mantenerle
			if( !IO.getKEEP() ){
				for (int maskIdx=1; maskIdx<=numberOfMasks;maskIdx++){
					boolean parent = true;
					IO.deleteFolder(outBench + "/mask" + maskIdx, parent);
					IO.deleteFolder(outLogo + "/mask" + maskIdx, parent);
				}
			}


		}

		long stopTime = System.currentTimeMillis();
		long elapsedTimeGlobal = stopTime - startTimeGlobal;
		System.out.println("full run time: " + ((int)elapsedTimeGlobal/10)/100.0 + " s");
		System.out.println("risultati contenuti in " + outFolder);
		if(WARNING_OOBIN)		Warning.OutOfBin();

	}//fine main

}
