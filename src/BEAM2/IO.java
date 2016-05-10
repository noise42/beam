package BEAM2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Random;

import org.apache.commons.io.FileDeleteStrategy;
import org.apache.commons.io.FileUtils;

public class IO {
	private static boolean KEEP=false;
	private static boolean INIT_SUMMARY=false;
	static void readInput(String path, ArrayList<Motif> inputSequences){
		FileReader f = null;
		try{
			//apro file
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader reader = new BufferedReader(f);
			String s,name="",seq="", firstName="", nuc="", dotB="";
			boolean start=true;
			boolean subopt = false;

			while ((s=reader.readLine())!=null){

				//System.out.println(s);
				if(start){
					//se non trova che all'inizio c'e' un > allora problema
					if((s.substring(0,1)).equals(">")){
						name=s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]","");
						start=false;
						//						if(firstSequence){ //se e' la prima subopt segnati il nome
						//							firstName = s.substring(1).split("�")[0];
						//							firstSequence = false;
						//						}

					}else{
						System.out.println("Check file format !\ns.substring(0,1) not equal '>' !!\n" + s.substring(0,1));
						System.exit(-1);
					}

				}else{ //non stiamo all'inizio
					String firstChar=s.substring(0,1);
					if(firstChar.equals(">")){ //se e' un intestazione fasta -> aggiungi la precedente
						//						System.err.println(nuc + "\n" + seq + "\n" + dotB);

						if(subopt){ //se nel passo precedente e' stata riconosciuta una subopt allora non aggiungere un nuovo
							//Motif ad inputSequences, ma riempi l'arrayList di String dell'ultimo messo
							inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);

						}else{
							//aggiungi la subopt -0- a inputSequences
							inputSequences.add(new Motif(name, nuc, dotB, seq,0,0));// (aggiunge quella precedente!!)
							firstName = name;
							inputSequences.get(inputSequences.size()-1).addSuboptMask();


						}



						name = s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]", "");
						seq="";
						nuc="";
						dotB="";


						if(name.equals(firstName)){//se la nuova sequenza e' una subopt di quella sopra
							subopt = true;
						}else{
							subopt = false;
						}

					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[acgurymkswbdhvntACGURYMKSWBDHVNT]*$")){ 
						//controllo per presenza di riga di sequenza
						nuc += s.trim().replaceAll("[\r\n]", "");
					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[.()]*$")){
						//controllo per presenza riga dotBracket
						dotB+= ((s.trim()).replaceAll("[\r\n]",""));
					}else{
						seq+= ((s.trim()).replaceAll("[\r\n]",""));
					}

				}
			}

			//quando esci segnati l'ULTIMA sequenza del file
			if(subopt){
				inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);
			}else{
				inputSequences.add(new Motif(name, nuc,dotB, seq,0,0));
				inputSequences.get(inputSequences.size()-1).addSuboptMask();

			}
		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	static public void scream(int status){
		System.err.println("AAAAAAARGH!");
		System.exit(status);
	}

	public static void deletePrevious(String myDirectoryPath) throws IOException{


		File dir = new File(myDirectoryPath);
		File[] directoryListing = dir.listFiles();
		if (directoryListing != null) {
			/*			for (File child : directoryListing) {
				// Do something with child
				deleteFolder(child);
			}
			 */
			//deleteFolder(dir);
		}
	}

	public static void deleteFolder(String myDirectoryPath) throws IOException{
		File fin = new File(myDirectoryPath);
		//deletes all child of folder, not the parent
		for (File file : fin.listFiles()) {
			FileDeleteStrategy.FORCE.delete(file);
		}   
	}

	public static void deleteFolder(String myDirectoryPath, boolean parent) throws IOException{
		//deletes folder
		File fin = new File(myDirectoryPath);
		File[] directoryListing = fin.listFiles();
		if(directoryListing != null){ 
			for (File file : fin.listFiles()) {
				FileDeleteStrategy.FORCE.delete(file);
			}   
			if (parent==true){
				FileDeleteStrategy.FORCE.delete(fin);
			}
		}


	}

	public static void moveBestRuns(String baseName, int mask, int bestRunIdx,
			String outDad) throws IOException {
		//TODO copy best benchmark in final folder, use bestRunIdx to manage filenames basename_m<mask>_run<bestRunIdx>
		String sourcename= outDad +"/mask"+ mask + "/" + baseName + "_m" + mask + "_run" + bestRunIdx + ".txt"; 
		String destname=outDad +"/motifs/" + baseName + "_m" + mask + "_run" + bestRunIdx + ".txt"; 

		File source = new File(sourcename);
		File dest = new File(destname);
		copyUsingFileChannels(source,dest);

	}

	public static void moveBestRunsLogo(String baseName, int mask, int bestRunIdx,
			String outDad) throws IOException {
		//TODO copy best benchmark in final folder, use bestRunIdx to manage filenames basename_m<mask>_run<bestRunIdx>
		String sourcename= outDad +"/mask"+ mask + "/" + baseName + "_m" + mask + "_run" + bestRunIdx + "_wl.fa"; 
		String destname=outDad +"/motifs/" + baseName + "_m" + mask + "_run" + bestRunIdx + "_wl.fa"; 

		File source = new File(sourcename);
		File dest = new File(destname);
		copyUsingFileChannels(source,dest);

	}

	private static void copyUsingFileChannels(File source, File dest) throws IOException {
		try {
			FileUtils.copyFile( source, dest);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void maskedExit() {
		// TODO Auto-generated method stub
		System.out.println("No unmasked sequences left");
		System.exit(0);
	}

	static void readInputForSearch(String path, ArrayList<Motif> inputSequences){
		//riempie inputSequences di ci� che c'� in path. 
		//DOVREBBE prendere anche inputSequences parzialmente riempiti
		FileReader f = null;
		try{
			//apro file
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader reader = new BufferedReader(f);
			String s,name="",seq="", firstName="", nuc="", dotB="";
			boolean start=true;
			boolean subopt = false;

			while ((s=reader.readLine())!=null){

				//System.out.println(s);
				if(start){
					//se non trova che all'inizio c'e' un > allora problema
					if((s.substring(0,1)).equals(">")){
						name=s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]","");
						start=false;
						//						if(firstSequence){ //se e' la prima subopt segnati il nome
						//							firstName = s.substring(1).split("�")[0];
						//							firstSequence = false;
						//						}

					}else{
						System.out.println("Check file format !\ns.substring(0,1) not equal '>' !!\n" + s.substring(0,1));
						System.exit(-1);
					}

				}else{ //non stiamo all'inizio
					String firstChar=s.substring(0,1);
					if(firstChar.equals(">")){ //se e' un intestazione fasta -> aggiungi la precedente
						//						System.err.println(nuc + "\n" + seq + "\n" + dotB);

						if(subopt){ //se nel passo precedente e' stata riconosciuta una subopt allora non aggiungere un nuovo
							//Motif ad inputSequences, ma riempi l'arrayList di String dell'ultimo messo
							inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);

						}else{
							//aggiungi la subopt -0- a inputSequences
							inputSequences.add(new Motif(name, nuc, dotB, seq,0,0));// (aggiunge quella precedente!!)
							firstName = name;
							inputSequences.get(inputSequences.size()-1).addSuboptMask();


						}



						name = s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]", "");
						seq="";
						nuc="";
						dotB="";


						if(name.equals(firstName)){//se la nuova sequenza e' una subopt di quella sopra
							subopt = true;
						}else{
							subopt = false;
						}

					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[acgurymkswbdhvntACGURYMKSWBDHVNT]*$")){ 
						//controllo per presenza di riga di sequenza
						//nuc += s.trim().replaceAll("[\r\n]", "");
					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[.()]*$")){
						//controllo per presenza riga dotBracket
						//dotB+= ((s.trim()).replaceAll("[\r\n]",""));
					}else{
						seq+= ((s.trim()).replaceAll("[\r\n]",""));
					}

				}
			}

			//quando esci segnati l'ULTIMA sequenza del file
			if(subopt){
				inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);
			}else{
				inputSequences.add(new Motif(name, nuc,dotB, seq,0,0));
				inputSequences.get(inputSequences.size()-1).addSuboptMask();

			}
		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	static void readInputForSearch(String path, ArrayList<Motif> inputSequences, boolean bgFlag){
		//riempie inputSequences di ci� che c'� in path. 
		//prende anche inputSequences parzialmente riempiti
		//VERSIONE background, appende _bg a tutte le sequenze, in modo da riconoscerle
		//per l'MCC (o altre analisi a posteriori)
		FileReader f = null;
		try{
			//apro file
			f = new FileReader(path);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader reader = new BufferedReader(f);
			String s,name="",seq="", firstName="", nuc="", dotB="";
			boolean start=true;
			boolean subopt = false;

			while ((s=reader.readLine())!=null){

				//System.out.println(s);
				if(start){
					//se non trova che all'inizio c'e' un > allora problema
					if((s.substring(0,1)).equals(">")){
						name=s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]","") + "_bg";
						start=false;
						//						if(firstSequence){ //se e' la prima subopt segnati il nome
						//							firstName = s.substring(1).split("�")[0];
						//							firstSequence = false;
						//						}

					}else{
						System.out.println("Check file format !\ns.substring(0,1) not equal '>' !!\n" + s.substring(0,1));
						System.exit(-1);
					}

				}else{ //non stiamo all'inizio
					String firstChar=s.substring(0,1);
					if(firstChar.equals(">")){ //se e' un intestazione fasta -> aggiungi la precedente
						//						System.err.println(nuc + "\n" + seq + "\n" + dotB);

						if(subopt){ //se nel passo precedente e' stata riconosciuta una subopt allora non aggiungere un nuovo
							//Motif ad inputSequences, ma riempi l'arrayList di String dell'ultimo messo
							inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);

						}else{
							//aggiungi la subopt -0- a inputSequences
							inputSequences.add(new Motif(name, nuc, dotB, seq,0,0));// (aggiunge quella precedente!!)
							firstName = name;
							inputSequences.get(inputSequences.size()-1).addSuboptMask();


						}



						name = s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]", "") + "_bg";
						seq="";
						nuc="";
						dotB="";


						if(name.equals(firstName)){//se la nuova sequenza e' una subopt di quella sopra
							subopt = true;
						}else{
							subopt = false;
						}

					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[acgurymkswbdhvntACGURYMKSWBDHVNT]*$")){ 
						//controllo per presenza di riga di sequenza
						//nuc += s.trim().replaceAll("[\r\n]", "");
					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[.()]*$")){
						//controllo per presenza riga dotBracket
						//dotB+= ((s.trim()).replaceAll("[\r\n]",""));
					}else{
						seq+= ((s.trim()).replaceAll("[\r\n]",""));
					}

				}
			}

			//quando esci segnati l'ULTIMA sequenza del file
			if(subopt){
				inputSequences.get(inputSequences.size()-1).getSequenceList().add(seq);
			}else{
				inputSequences.add(new Motif(name, nuc,dotB, seq,0,0));
				inputSequences.get(inputSequences.size()-1).addSuboptMask();

			}
		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}



	public static void keepRuns() {
		//metodo attivato se c'� la flag -k
		IO.KEEP=true;
	}
	public static boolean getKEEP(){
		return IO.KEEP;
	}

	public static void writeGaussian(String pathBench, double gaussMean, double gaussVar) {
		Writer writer = null;
		try {

			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(pathBench), "utf-8"));


			writer.write("mean\t" +gaussMean);
			writer.write("\n");
			writer.write("var\t" + gaussVar);

		} catch (IOException ex) {
		} finally {
			try {writer.close();} catch (Exception ex) {}
		}

	}
	
	public static void writeBenchmark(String pathBench, MotifHandler mh, boolean escaped, int contatore, int minSteps, int passiTot){
		//write Benchmark dato il nome del file, l mh, e i parametri di uscita (sostituire anche in run) TODO
		Writer writer = null;
		try {

			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(pathBench), "utf-8"));


			writer.write(mh.toString());
			writer.write("\n#otherMatches\n");
			writer.write(mh.otherMatches());

			writer.write("\n#Seq PSSM\n");
			writer.write(mh.toStringSeq());			



			writer.write("\n#PSSM\n");
			writer.write(mh.printPSSM());
			writer.write("\n#score\t" + mh.getScore());
			writer.write("\n#seq\t" + mh.cardinality());
			writer.write("\n#width\t" + mh.getMotifWidth());
			//					if (mh.getObjectMotif(0).getNucleotides().equals("") == false){
			//						writer.write("\n#nucCons\t" + PSSM.checkBonusSeq(bestMh));
			//					}else{
			//						writer.write("\n#nucCons\tNULL");
			//					}
			writer.write("\n#escape\t" + escaped);
			writer.write("\n#onStep\t" + contatore);
			writer.write("\n#minSteps\t" + minSteps);
			writer.write("\n#maxSteps\t" + passiTot);


		} catch (IOException ex) {
		} finally {
			try {writer.close();} catch (Exception ex) {}
		}

	}

	public static void writeFailedSearchBenchmark(String pathBench) {
		//write null Benchmark after search if no Negatived were put in
		//In realt� una volta che ci si aggiunge l'input al bg non c'� necessit� di sta cosa
		//perch� ci sar� sempre almeno una sequenza (quelle che compongono il mh)
		Writer writer = null;
		String MESSAGE = "No sequences were found. This means problems. Contact the coders @ marco.pietrosanto at uniroma2.it";
		try {

			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(pathBench), "utf-8"));


			writer.write(MESSAGE);

		} catch (IOException ex) {
		} finally {
			try {writer.close();} catch (Exception ex) {}
		}

	}

	public static SearchObj readSearchOutput(String searchFilePath,int totalInput,int totalBackground) {
		//legge formato .search.txt, ne ricava dati utilizzati da
		
		FileReader f = null;
		String s="";
		//farsi tutta la PSSM mi pare esagerato, va bene una matrice di PSSMcell
		ArrayList<ArrayList<PSSMCell>> PSSMdata = new ArrayList<ArrayList<PSSMCell>>();
		//inizializzo dati per l'mcc
		int TP,FP,TN,FN;
		TP=FP=TN=FN=0;
		double mediaPart=0.0;
		try{
			//apro file
			f = new FileReader(searchFilePath);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}

		try{
			BufferedReader reader = new BufferedReader(f);
			//riconoscere flag messa al bg
			boolean alignmentZone=true;
			boolean PSSMzone=false;
			while ((s=reader.readLine())!=null){
				if (s.trim().equals("")){
					//questa � generale
					PSSMzone=false;
					alignmentZone=false;
				}
				//le operazioni vanno prima delle flag, altrimenti leggo anche #PSSM etc.
				if (alignmentZone){
					//segna se � TP o FP
					if (s.split("\t")[1].split("[$]")[0].endsWith("_bg")){
						FP++;
					}else{
						TP++;
						//calcolo media parziali
						mediaPart += Double.parseDouble(s.split("\t")[5]);
					}
				}
				if (PSSMzone){
					//riempi PSSMdata
					ArrayList<PSSMCell> tmpAL = new ArrayList<PSSMCell>();
					for (String chunk: s.trim().split("\t")){
						PSSMCell c = new PSSMCell(chunk.charAt(0), Double.parseDouble(chunk.split(":")[chunk.split(":").length-1]));
						tmpAL.add(c);
					}
					PSSMdata.add(tmpAL);
				}

				//flags in mezzo al ciclo
				if (!alignmentZone && s.startsWith("#PSSM")){
					PSSMzone=true;
				}

				//fine flags

			}
			FN = totalInput - TP;
			TN = totalBackground - FP;
			mediaPart /= TP;
			
		}catch(IOException ioException){ioException.printStackTrace();}

		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		SearchObj searchObj = new SearchObj(PSSMdata, TP, FP, TN, FN, mediaPart);
		return searchObj;
	}


	public static void appendToSummary(String summaryName, SearchObj sObj, SummarizedData sumData, long elapsedTime){
		
		//DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.ITALIAN);
		//otherSymbols.setDecimalSeparator('.');
		//NumberFormat formatter = new DecimalFormat("0.##E0"); 

		try {
			FileWriter writer = new FileWriter(summaryName, true);

			if(!INIT_SUMMARY){
				writer.write("#input=" + sumData.getInputFile());
				writer.write("\n");
				writer.write("#input_size=" + sumData.getInputSize());
				writer.write("\n");
				writer.write("#background=" + sumData.getBgFile());
				writer.write("\n");
				writer.write("#background_size=" + sumData.getBgSize());
				writer.write("\n");
				writer.write("#meanLength=" + sumData.getMl());
				writer.write("\n");
				writer.write("#meanStructureContent=" + sumData.getMsc());
				writer.write("\n");
				writer.write("#bin=" + sumData.computeBin());
				writer.write("\n");
				INIT_SUMMARY=true;
			}
			writer.write("\n");
			writer.write("#motif=" + sumData.getMask());
			writer.write("\n");
			writer.write(sumData.getConsensus());
			writer.write("#qBEAR\n");
			writer.write(sumData.getqBEARConsensus());


//			writer.write("MCC=" + sumData.getMCC());
//			writer.write("\n");
			writer.write("pvalue=" + sObj.getPvalue());
			writer.write("\n");			
			writer.write("score=" + sumData.getScore());
			writer.write("\n");
			writer.write("coverage=" + sumData.getCoverage());
			writer.write("\n");
			writer.write("fall-out=" + sumData.getFallout());
			writer.write("\n");
			//writer.write("accuracy=" + sumData.getAcc());
			//writer.write("\n");
//			writer.write("TP=" + sObj.getTP());
//			writer.write("\n");
//			writer.write("FP=" + sObj.getFP());
//			writer.write("\n");
//			writer.write("TN=" + sObj.getTN());
//			writer.write("\n");
//			writer.write("FN=" + sObj.getFN());
//			writer.write("\n");
			writer.write("runtime(s)="+ (((int)elapsedTime/10)/100.0));
			writer.write("\n");




			writer.close();
		} catch (IOException ex) {
			System.err.println("IOException: " + ex.getMessage());
		}
	}

	public static void appendToBenchmark(String baseName, int mask, int bestRunIdx) {
		// TODO appende i dati della best run (consensus, score, mcc) al file di benchmark

	}

	public static String defaultBgCreation(String bin, int size, String outFolder) {
		//Crea file di background pescando da bin, max 250 (come da RFtest) 
		//restituisce nome del file creato, da passare a SearchAfterRun
		String BGSRC="/defaultBg/";
		String bgFile=BGSRC+ "noise_rfam_list" + bin + ".fa";
		InputStream stream = IO.class.getResourceAsStream(bgFile);
		
		//file di output
		String outBgFileName=outFolder + "/autoBG.fa";
		//arrayList in cui saranno storati gli rna del bin
		ArrayList<Motif> inputSequences = new ArrayList<Motif>();

		if (size>250) size=250; //massima eterogeneit� consentita da Rfam
		
		FileReader f = null;
/*		try{
			//apro file
			//f = new FileReader(bgFile);
		}catch (IOException ioException){
			ioException.printStackTrace();
		}
*/
		try{
			BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
			String s,name="",seq="", firstName="", nuc="", dotB="";
			boolean start=true;

			while ((s=reader.readLine())!=null){

				//System.out.println(s);
				if(start){
					//se non trova che all'inizio c'e' un > allora problema
					if((s.substring(0,1)).equals(">")){
						name=s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]","");
						start=false;
						//						if(firstSequence){ //se e' la prima subopt segnati il nome
						//							firstName = s.substring(1).split("�")[0];
						//							firstSequence = false;
						//						}

					}else{
						System.out.println("Check file format !\ns.substring(0,1) not equal '>' !!\n" + s.substring(0,1));
						System.exit(-1);
					}

				}else{ //non stiamo all'inizio
					String firstChar=s.substring(0,1);
					if(firstChar.equals(">")){ //se e' un intestazione fasta -> aggiungi la precedente
						//						System.err.println(nuc + "\n" + seq + "\n" + dotB);

						
						//aggiungi la subopt -0- a inputSequences
						inputSequences.add(new Motif(name, nuc, dotB, seq,0,0));// (aggiunge quella precedente!!)
						firstName = name;
						inputSequences.get(inputSequences.size()-1).addSuboptMask();


					



						name = s.substring(1).split("\\$")[0].split("pos")[0].replaceAll("[\r\n]", "");
						seq="";
						nuc="";
						dotB="";




					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[acgurymkswbdhvntACGURYMKSWBDHVNT]*$")){ 
						//controllo per presenza di riga di sequenza
						nuc += s.trim().replaceAll("[\r\n]", "");
					}else if(s.replaceAll("[\r\n]", "").trim().matches("^[.()]*$")){
						//controllo per presenza riga dotBracket
						dotB+= ((s.trim()).replaceAll("[\r\n]",""));
					}else{
						seq+= ((s.trim()).replaceAll("[\r\n]",""));
					}

				}
			}

			//quando esci segnati l'ULTIMA sequenza del file

			inputSequences.add(new Motif(name, nuc,dotB, seq,0,0));
			inputSequences.get(inputSequences.size()-1).addSuboptMask();


		}catch(IOException ioException){ioException.printStackTrace();}

		/*
		try {
			f.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		*/
		
		//A questo punto inputSequences contiene tutto il Bg, 
		//ma vanno selezionati i motif presi a random
		
		ArrayList<Motif> bgSubset = new ArrayList<Motif>();
		
		generateRandomSubset(inputSequences, bgSubset, size);
		
		//scrivo il file finale in outBgFileName
		
		try {
			FileWriter writer = new FileWriter(outBgFileName, true);

			for(Motif m: bgSubset){
				writer.write(">" + m.getName());
				writer.write("\n");
				writer.write(m.getNucleotides());
				writer.write("\n");
				writer.write(m.getDotBracket());
				writer.write("\n");
				writer.write(m.getSequence());
				writer.write("\n");
			}

			writer.close();
		} catch (IOException ex) {
			System.err.println("IOException: " + ex.getMessage());
		}
		
		return outBgFileName;
	}

	private static void generateRandomSubset(ArrayList<Motif> inputSequences,
			ArrayList<Motif> bgSubset, int size) {
			//Pesca un numero di sequenze pari a size, rispettando le diverse famiglie
			//(...baciamo le mani)
		
		Collections.shuffle(inputSequences, new Random(System.currentTimeMillis()));
		for(int i=0;i<size;i++){
			bgSubset.add(inputSequences.get(i));
		}
		
	}


}
