package BEAM2;

//MULTI RUN BRANCH
//___
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;

public class BEAM2 {
	private static boolean acceptScore(double scoreVecchio, double scoreNuovo, double temperature) {
		if (scoreNuovo-scoreVecchio >0) System.err.println("SPOTTED SOME DISCREPANCIES!! Score mandato in"
				+ "accettazione quando non doveva");
		double prob = Math.exp(1*(scoreNuovo - scoreVecchio)/temperature);
		if (temperature == 0){
			prob = 0;
		}
		if(Math.random() < prob ){
			return true;
		}
		else{
			return false;
		}
	}

	/**
	 * @param args
	 */
	public static MotifMemory run(String[] args, int idx, int mask, String output, ArrayList<Motif> inputSequences) {
		long startTime = System.currentTimeMillis();

		//ArrayList<Motif> inputSequences=new ArrayList<Motif>();

		MotifUtilities.reshapeInput(inputSequences);
		
		int escape = 0;

		double[][]mbr=MBR.getmbr();	

		MotifManager manager=new MotifManager();

		//default
		int motifWidth=10;
		int startingNoOfSequences = 3;
		int model_limit = 100;

		double temperatureStart=100;
		double coolingRate=0.001;//These two variables can be fixed values or can be taken as input
		double minTemp = 0.001;


		int widthUpperBound = 100;

		int minSteps = 10000;
		double passiTot=15000; 

		int escapeCondition = 1;



		final int N = 5; //numero medio di operazioni per sequenza:controlla il numero minimo di passi
		//CommandLine Parser 
		//TODO pure tutti questi comandi sarebbe meglio che li leggesse una volta sola
		CommandLineParser clp = new CommandLineParser(args);
		//initialized CLP, the hash is created, check with - GetValue(key), ContainsKey(key)

		//do not touch
		double maximumScoreEver=-9999;
		boolean escaped = false;
		double temperature = temperatureStart;
		boolean clean = true; //parziali > 50% media
		boolean superClean = true; //parziali >0
		boolean ultraClean = true; //parziali >0

		boolean branching = true;
		boolean weHaveSubopt = false;

		//-----LETTURA COMMAND LINE ARGUMENTS (CLA)
		if(clp.containsKey("s")) minSteps=Integer.parseInt(clp.getValue("s")); //numero min di passi prima di fermarsi
		if(clp.containsKey("W")) widthUpperBound=Integer.parseInt(clp.getValue("W")); //limite superiore larghezza motif
		if(clp.containsKey("w")) motifWidth=Integer.parseInt(clp.getValue("w")); //starting motif width
		if(clp.containsKey("T")) temperature=Double.parseDouble(clp.getValue("T")); //starting T
		if(clp.containsKey("r")) coolingRate=Double.parseDouble(clp.getValue("r")); //cooling rate, precision

		if(clp.containsKey("n")) model_limit=Integer.parseInt(clp.getValue("n")); //sampling iniziale
		if(clp.containsKey("o")) output=clp.getValue("T"); //output folder
		if(clp.containsKey("C") && clp.getValue("C").equals("1")) {clean=true;superClean=false;ultraClean=false;}
		else if(clp.containsKey("C") && clp.getValue("C").equals("2")){clean=true;superClean=true;ultraClean=false;}
		else if(clp.containsKey("C") && clp.getValue("C").equals("3")){clean=true;ultraClean=true;}
		if(clp.containsKey("b")) {
			if (clp.getValue("b").equals("F")) branching=false;
		}
		if(clp.containsKey("m") && clp.getValue("m").equals("T")) weHaveSubopt = true;
		//----------
		//		percentage =startingNoOfSequences/(double)inputSequences.size();
		
		if (!branching){
			MotifUtilities.setBranchingZero(mbr);
			if (Debug.VERBOSE){
				System.out.println("Setting :-: equal to 0...mbr[:][:]: " + BEARManager.substitutionScore(':', ':', mbr));
			}
		} else {
			if (Debug.VERBOSE){
				System.out.println("branch to branch substitution considered: " + BEARManager.substitutionScore(':', ':', mbr));
			}
		}


		coolingRate = (temperatureStart - minTemp)/minSteps; //in pratica quando la T va a 0.0


		//-----QUI FINISCONO I PARAMETRI DI INPUT-----

		//		System.out.println(motifWidth + " " + inputSequences.size() + " " + percentage);
		MotifHandler mh = null;
		MotifHandler bestMh = new MotifHandler();
		int mW=0;

		if ( mask == 1 && idx == 1){
			//System.out.println(inputSequences.size());
			mh=manager.initialise(motifWidth,inputSequences,mbr, startingNoOfSequences);//genera la prima soluzione casuale		
		}else{
			mh=manager.initialise2(motifWidth,inputSequences,mbr,startingNoOfSequences);
		}
		if (Debug.VERBOSE) {System.out.println("---GENERATED STARTING MH---\n"+mh.toString());}
		//---------------------
		PSSM.computePSSM(mh, inputSequences);
		mh.setMotifWidth();		

		if (Debug.VERBOSE){
			if (mh.getObjectMotif(0).getNucleotides().equals("")){
				System.err.println("Nucleotide sequences not found");
			}
		}
		
		if (Debug.ON) Debug.checkRNA(mh);
		MotifUtilities.computePartials(mh.getPSSM(), mh, mbr, mh.getMotifWidth());
		MotifHandler.computeScore(mh);//genera primo punteggio
		double prevScore = mh.getScore();
		int contatore=0;

		//per avere almeno N operazioni su ogni sequenza in media (7 operazioni, X subopt)
		//CLA
		if(clp.containsKey("S")) {
			passiTot=Double.parseDouble(clp.getValue("S")); //force passiTot
		}else{


			//-----------
			passiTot = inputSequences.size()*7*mh.getObjectMotif(0).getSequenceList().size()*N; 		

			if (passiTot < minSteps){
				passiTot = minSteps + inputSequences.size()*7*mh.getObjectMotif(0).getSequenceList().size();
			}			
			if (passiTot > 15000){
				passiTot=15000;
			}
		}
		

		escapeCondition = Math.max((int)(passiTot*.005),1000); //anche in un dataset piccolo richiediamo almeno 1000 steps per scappare



		//se dopo aver in media agito con ogni perturbazione su tutte le sequenze una volta non cambia nulla, allora escape, ma solo se sono stati fatti almeno minSteps

		if (Debug.VERBOSE) {
			System.out.println("..............STARTING MOTIF..............\n");
			System.out.println(mh.toString());
			System.out.println("Score : " + mh.getScore());
			System.out.println("Length : " + mh.getMotifWidth());
			System.out.println("#Sequences : " + mh.cardinality() +"\n");
		}

		/*---------*/
		double[][] media;
		media = new double[][]{{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
		int[][] count;
		count = new int[][]{{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
		int[][] acceptSecondCount;
		int[][] rejectCount;
		acceptSecondCount = new int[][]{{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};
		rejectCount = new int[][]{{0,0,0,0,0,0,0},{0,0,0,0,0,0,0},{0,0,0,0,0,0,0}};

		/*---------*/





		//

		double currentScore=mh.getScore();
		
		//-----INIZIO CICLO -----

		while((escape < escapeCondition || contatore < minSteps) && contatore < passiTot ){
			//se anche e' stata raggiunta la condizione di escape non fermarti prima di minSteps, in ogni caso esci dopo passiTot
			if (Debug.VERBOSE) {System.out.println("\n>>>>>>>>>>>>>>>ITERATION : "+contatore+ " di " + passiTot + "<<<<<<<<<<<<<<<");}
			currentScore=mh.getScore();
			//System.out.println(perturbatedMotif.motifLength + "\t" + m.motifLength);

			if(mh.cardinality() >= 1){
				mh.setMotifWidth();
			}else{
				//				System.err.println(mh.cardinality());
			}

			//-----PERTURBO IL MOTIVO, TENGO LA SOLUZIONE PRECEDENTE-----
			int operation = manager.perturbateMotif(mh,inputSequences,mbr, widthUpperBound, weHaveSubopt, model_limit);

			//-----NEL CASO NORMALE CALCOLA PSSM E SCORE-----
			if(mh.cardinality()<2){
				//-----SE IL MH CONTIENE 1 SEQUENZE INVECE FAI UN ADD-----		
				currentScore=0;
				mh.setMotifWidth(motifWidth);
				operation = manager.perturbateMotif(mh,inputSequences,mbr, widthUpperBound, weHaveSubopt, model_limit);
			}
			if (operation == -1){
				break;
			}
			//SE CI SONO PROBLEMI DECOMMENTARE QUESTO SETMOTIFWIDTH ALERT ALERT ALERT!
			//mh.setMotifWidth();

			if (Debug.VERBOSE){
				System.out.println(mh.toString());
			}
			PSSM.computePSSM(mh, inputSequences);
			MotifUtilities.computePartials(mh.getPSSM(), mh, mbr, mh.getMotifWidth());
			MotifHandler.computeScore(mh);

			/*Test per peso singole operazioni*/
			double delta = mh.getScore() - currentScore;

			if (contatore < minSteps){
				media[(int)(contatore/(minSteps/3.0))][operation] += delta;
				count[(int)(contatore/(minSteps/3.0))][operation] += 1;
			}
			/*-------------------------------*/
			double newScore = mh.getScore();

			if(newScore>=currentScore){
				if (Debug.VERBOSE) {System.out.println("---ACCEPTED---\n");}
				currentScore=newScore;

			}else if(acceptScore(currentScore,newScore, temperature)){
				if (Debug.VERBOSE) {System.out.println("---ACCEPTED SECOND CHANCE---\n");}				
				currentScore=mh.getScore();

				acceptSecondCount[(int)(contatore/(passiTot/3))][operation] += 1;
			}else{
				if (Debug.VERBOSE) {System.out.println("\n---REJECTED---\n");}
				MotifManager.ctrlZ(operation, mh, inputSequences);//il ctrlZ non riprende il vecchio score!!! ora si
				mh.getPSSM().setScore(currentScore); //aggiunto : maxScore e' semplicemente lo score del giro precedente

				if (contatore<minSteps-1){
					rejectCount[(int)(contatore/minSteps*3)][operation] += 1;
				}

			}


			if (Debug.VERBOSE) {	System.out.println("---cooling down---\n");
			}
			if(temperature > minTemp){
				temperature-=coolingRate;
			}else{
				temperature = 0;
			}


			//----------------------------TEST MAX SCORE----------------------------


			if(newScore >= prevScore){
				//se migliora lo score, non aumentare l'escape, anzi azzeralo. A meno che non sia aumentato di troppo poco.
				maximumScoreEver = newScore;

				bestMh.MotifHandlerClone(mh);
				bestMh.getPSSM().setScore(maximumScoreEver); 


				int threshold = 1;
				if(Math.abs(mh.getScore() - prevScore) < threshold){ //TODO da stimare meglio
					escape++;
					//System.err.println("escape in: " + (escapeCondition - escape) + " passi");
				}else{
					escape = 0;
				}
			}
			prevScore = newScore;
			//----------------------------TEST MAX SCORE END----------------------------


			if(mh.cardinality() >= 1){
				mh.setMotifWidth();
			}
			
			contatore++;

			if (Debug.VERBOSE) {
				System.out.println("Temperature : " + temperature);
				System.out.println("Score : " + mh.getScore());
				System.out.println("Length : " + mh.getMotifWidth());
				System.out.println("#Sequences : " + mh.cardinality() +"\n");
			}else if(contatore%(int)(passiTot/10) == 0 ){
				System.out.println(Math.round(contatore/passiTot*100) + "%");
			}
			mW=bestMh.getMotifWidth();
		}
		System.out.println("\n-------------END OF COMPUTATION-------------\n");

		bestMh.adjustEndIndexes(mW);
		if(escape >= escapeCondition) escaped = true;

		
		if (bestMh.cardinality() == 0){
			bestMh.MotifHandlerClone(mh);
			PSSM.computePSSM(bestMh, inputSequences);
			bestMh.getPSSM().setScore(mh.getScore());
		}else{
			PSSM.computePSSM(bestMh, inputSequences);
		}
		
		if (Debug.VERBOSE){
		
			System.out.println("BestMh PSSM --------------------");
			System.out.println(bestMh.printPSSM());

		}


		MotifUtilities.computePartials(bestMh.getPSSM(), bestMh, mbr, bestMh.getMotifWidth());

		//Pulizia --------- (--clean <- non funge? -C? non funge?)
		System.out.println("clean: " + clean);

		if (clean){
			System.out.println("\n[-------------CLEANING INITIATED-------------]\n");

			double thr = 0.0;
			if(superClean){
				thr = 0.5;
			}else if(ultraClean){
				thr= .99;
			}
			MotifUtilities.clean(bestMh, thr);

			PSSM.computePSSM(bestMh, inputSequences); //rifaccio dopo il clean
			motifWidth=bestMh.getMotifWidth();
			MotifUtilities.computePartials(bestMh.getPSSM(), bestMh, mbr, motifWidth);

			//System.out.println("Motif Score : "+mh.getScore());
			System.out.println("Motif Score : "+bestMh.getScore());

			MotifHandler.computeScore(bestMh);

			//System.out.println("Motif Score after Cleaning : "+mh.getScore());
			System.out.println("Motif Score after Cleaning : "+bestMh.getScore());
			System.out.println("\n[-------------CLEANING COMPLETED-------------]\n");
		}


		if(mh.cardinality() != 0){
			mh.setMotifWidth();}


		//Media Delta Score
		if (Debug.ON){ 
			System.out.println("\n-------------Medie Delta Score-------------\n");
			for(int parte = 0; parte < 3; parte++){
				for(int k = 0; k < media[parte].length; k++){
					media[parte][k] /= count[parte][k];
				}
			}	

			for (int phase = 0; phase<=2; phase++){
				if(phase == 0){
					System.out.println("\nfase high T-----");
				}else if(phase == 1){
					System.out.println("\nfase mid T-----");
				}else{
					System.out.println("\nfase low T-----");
				}
				for (int op = 0; op<=5; op++){
					if (op == 0){
						System.out.print("media add: ");
					}else if(op == 1){
						System.out.print("media remove: ");
					}else if(op == 2){
						System.out.print("media shift: ");
					}else if(op == 3){
						System.out.print("media recalculate: ");
					}else if(op == 4){
						System.out.print("media expand: ");
					}else if(op == 5){
						System.out.print("media shrink: ");
					}
					System.out.println(media[phase][op] + "\t#operazioni: " + count[phase][op] +
							"\taccept: " + (count[phase][op]-acceptSecondCount[phase][op]-rejectCount[phase][op]) + "\taccept2nd: " + acceptSecondCount[phase][op] + 
							"\treject: " + rejectCount[phase][op]);
				}
			}
		}
		
		System.out.println("\n----- final local alignment -----\n");
		//riportare comunque il maximum score ever
//		if(maximumScoreEver > mh.getScore()){
//			
//			System.out.println(bestMh.toString());
//			if (bestMh.getObjectMotif(0).getNucleotides().equals("") == false){
//				System.out.println(bestMh.toStringSeq());
//			}
//			System.out.println("\nSequences: " + bestMh.cardinality());
//			System.out.println("\nbestMh score: " + maximumScoreEver);
//
//		}else{
//			bestMh.getPSSM().setScore(mh.getScore());
//			System.out.println(bestMh.toString());
//			if (bestMh.getObjectMotif(0).getNucleotides().equals("") == false){
//				System.out.println(bestMh.toStringSeq());
//			}
//			System.out.println("\n\nn.Sequenze: " + bestMh.cardinality());
//			System.out.println("\n\nscore bestMh: " + bestMh.getScore());			
//		}

//		if (mh.getObjectMotif(0).getNucleotides().equals("") == false){
//			System.out.println("\n#nucleotide conservation (given as the fraction of conserved columns (>70%)\t" + PSSM.checkBonusSeq(bestMh));
//		}

		String fileBench="lastrun_m" + mask + "_run"  + idx ;
		if (clp.containsKey("f")){
			String[] tmp =  clp.getValue("f").split("/");
			fileBench = tmp[tmp.length-1].split("\\.")[0]+ "_m" + mask + "_run"  + idx; //
		}

		//WebLogo

		String webLogoOut = output + "/webLogoOut/mask"+ mask +"/" + fileBench + "_wl.fa"; //output da dare a weblogo
		boolean webLogo=false;
		if(webLogo){
			System.out.println(bestMh.toString(webLogo,true));
		}	

		Writer writer = null;

		try {
			writer = new BufferedWriter(new OutputStreamWriter(

					new FileOutputStream(webLogoOut), "utf-8"));


			writer.write(bestMh.toString(webLogo, true));

		} catch (IOException ex) {
		} finally {
			try{writer.close();} catch (Exception ex) {}
		}

		//write Benchmark

		String pathBench = output + "/benchmark/mask" + mask + "/" + fileBench + ".txt";
		try {

			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(pathBench), "utf-8"));


			writer.write(bestMh.toString());
			writer.write("\n#otherMatches\n");
			writer.write(bestMh.otherMatches());

			if (bestMh.getObjectMotif(0).getNucleotides().equals("") == false){
				writer.write("\n#Seq PSSM\n");
				writer.write(bestMh.toStringSeq());			
			}


			writer.write("\n#PSSM\n");
			writer.write(bestMh.printPSSM());
			writer.write("\n#score\t" + bestMh.getScore());
			writer.write("\n#seq\t" + bestMh.cardinality());
			writer.write("\n#width\t" + bestMh.getMotifWidth());
//			if (mh.getObjectMotif(0).getNucleotides().equals("") == false){
//				writer.write("\n#nucCons\t" + PSSM.checkBonusSeq(bestMh));
//			}else{
//				writer.write("\n#nucCons\tNULL");
//			}
			writer.write("\n#escape\t" + escaped);
			writer.write("\n#onStep\t" + contatore);
			writer.write("\n#minSteps\t" + minSteps);
			writer.write("\n#maxSteps\t" + passiTot);


		} catch (IOException ex) {
		} finally {
			try {writer.close();} catch (Exception ex) {}
		}



		//file da passare al jar che calcola lo score, per i test

		String pathScoreTest = output + "/scoreTest/mask" + mask + "/" + fileBench + ".txt";

		try {

			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(pathScoreTest), "utf-8"));


			writer.write(bestMh.toString2());


		} catch (IOException ex) {
		} finally {
			try {writer.close();} catch (Exception ex) {}
		}

		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.println("run time: " + ((int)elapsedTime/10)/100.0 + " s");

		return new MotifMemory(bestMh.getScore(), bestMh);

	}

	public static void applyMask(MotifHandler mh){
		for (Motif m: mh.getListMotif()){
			m.addMask(m.getMotifStart(), m.getMotifEnd());
			//inputSequences.add(m); //no perch� alla fine deve mascherare con la migliore ma riepire dell'ultima
		}
	}

	public static void refillInput(MotifHandler mh, ArrayList<Motif> inputSequences){
		if (Debug.ON){
			System.out.println("input size: "+ inputSequences.size());
			System.out.println("mh size: " + mh.cardinality());
			System.out.println();
		}
		for (Motif m: mh.getListMotif()){
			inputSequences.add(m);
		}
	}

	public static void main(String[] args) throws IOException {

		//inputSequences
		ArrayList<Motif> inputSequences=new ArrayList<Motif>();


		int numberOfRuns = 1;
		int numberOfMasks=1;
		String baseName = "lastrun";
		//CommandLine Parser 
		CommandLineParser clp = new CommandLineParser(args);
		//initialized CLP, the hash is created, check with - GetValue(key), ContainsKey(key)
		if(clp.containsKey("R")){
			numberOfRuns = Integer.parseInt(clp.getValue("R"));
		}
		if(clp.containsKey("f")){
			String[] tmp =  clp.getValue("f").split("/");
			baseName = tmp[tmp.length-1].split("\\.")[0]; //
			IO.readInput(clp.getValue("f"), inputSequences);
		}


		String outFolder="risultati/" + baseName;
		final File dir = new File(outFolder);
		if(dir.exists()){

			System.out.println("My scans tell me I already own an output with this name -"+baseName+ "-, want me to overwrite it?");

			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(System.in));
			String response = bufferedReader.readLine();
			//String response = console.readLine("y/n");
			if(response.toLowerCase().equals("y")){
				System.out.println("OVERWRITING...");
			}else{
				System.out.println("CLOSING...");
				System.exit(0);
			}
		}
		String outBench=outFolder+"/benchmark";
		String outScore=outFolder+"/scoreTest";
		String outLogo=outFolder+"/webLogoOut";
		dir.mkdirs();
		new File(outBench).mkdirs();
		new File(outScore).mkdirs();
		new File(outLogo).mkdirs();


		for(int mask=1; mask<=numberOfMasks; mask++){
			for(int i=1; i<=numberOfRuns; i++){
				run(args, i,mask, outFolder, inputSequences);
				System.out.println("run " + i + " completed");
			}
		}
	}


}

