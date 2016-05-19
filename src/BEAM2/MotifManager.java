package BEAM2;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class MotifManager {

	public MotifHandler initialise(int motifLength, ArrayList<Motif> inputSequences,double[][] m, int startingNo) {
		System.out.println("---INITIALIZING ---");
		MotifHandler tmp=new MotifHandler();
		Motif startingMotif;
		//shuffle the input sequences
		Collections.shuffle(inputSequences, new Random(System.currentTimeMillis()));
		//take the first in the stack and add it to the motifHandler tmp
		tmp.addMotif(inputSequences.get(0));
		startingMotif=tmp.getObjectMotif(0);
		//rimuovi la sequenza presa dallo stack di input
		inputSequences.remove(startingMotif);

		Random generator = new Random();
		//System.out.println(startingMotif.getSequence() + "\t" + motifLength);
		int roll = generator.nextInt(startingMotif.getSequence().length()-motifLength);//

		if(roll+motifLength > startingMotif.getSequence().length()) System.err.println("tentato motifStart maggiore del massimo consentito: " 
				+ roll + ", " + "lunghezza seq: " + startingMotif.getSequence().length());

		tmp.getObjectMotif(0).setMotifStart(roll, motifLength);
		tmp.getObjectMotif(0).setMotifEnd(roll+motifLength, motifLength);

		tmp.setMotifWidth();

		int shortestLength = 9999999;

		if(startingNo < 1){
			startingNo = 1;
		}
		for (int i = 0; i < (int)startingNo ; i++) {//da decidere
			Motif seq = new Motif();
			seq = inputSequences.get(0); //get Motif i
			int res=BEARManager.searchInitialMotif(startingMotif.extractMotifFromSequence(), seq ,m,  motifLength);

			if(res+motifLength > seq.getSequence().length()) System.err.println("tentato motifStart maggiore del massimo consentito: " 
					+ res + ", " + "lunghezza seq: " + seq.getSequence().length());

			seq.setMotifStart(res);
			seq.setMotifEnd(res+motifLength);
			tmp.addMotif(seq);
			inputSequences.remove(seq);
			int tmpShortest = seq.getSequence().length();
			if (tmpShortest < shortestLength){
				shortestLength = tmpShortest;
			}
		}
		tmp.setShortestSequenceLength(shortestLength);
		return tmp;
	}

	public MotifHandler initialise2(int motifLength,
			ArrayList<Motif> inputSequences, double[][] m, int startingNo) {
		//This method is used during the second and subsequent runs. It is pretty much the same
		//as initialise apart from it does not initialize masks
		System.out.println("---INITIALIZING __ ---");
		MotifHandler tmp=new MotifHandler();
		Motif startingMotif;


		Random generator = new Random();
		int roll = 99999;
		///*
		boolean goodSeq=false;
		boolean maskGoFlag=false;
		int maximumTries = inputSequences.size()*inputSequences.get(0).getSequenceList().size()*10;
		//*10 per aumentare le possibilit� di prenderne una che non e' stata presa
		int counter = 0;

		//shuffle the input sequences
		Collections.shuffle(inputSequences, new Random(System.currentTimeMillis()));
		//take the first in the stack and add it to the motifHandler tmp
		tmp.addMotif(inputSequences.get(0));
		startingMotif=tmp.getObjectMotif(0);

		while(!goodSeq){
			int rollCounter=0;
			while(!maskGoFlag){


				roll = generator.nextInt(startingMotif.getSequence().length()-motifLength);//
				if(roll+motifLength > startingMotif.getSequence().length() | roll < 0) {
					System.err.println("tentato motifStart maggiore del massimo consentito: " 
							+ roll + ", " + "lunghezza seq: " + startingMotif.getSequence().length());
				}

				for (int i=roll; i<=(roll+motifLength); i++){ //verifico tutte le posizioni coperte	
					maskGoFlag= !tmp.getObjectMotif(0).isMasked(i);
					if (!maskGoFlag) break;
				}
				if (maskGoFlag) break;

				rollCounter++;
				if (rollCounter > startingMotif.getSequence().length()){
					break;
				}

			}


			if (maskGoFlag) break;


			//else shuffle the input sequences
			Collections.shuffle(inputSequences, new Random(System.currentTimeMillis()));
			//take the first in the stack and add it to the motifHandler tmp
			tmp.removeMotif(0);
			tmp.addMotif(inputSequences.get(0));
			startingMotif=tmp.getObjectMotif(0);

			counter++; //se devi riniziare il ciclo aumenta il counter
			if (counter > maximumTries){
				IO.maskedExit();
			}

		}

		tmp.getObjectMotif(0).setMotifStart(roll, motifLength);
		tmp.getObjectMotif(0).setMotifEnd(roll+motifLength, motifLength);
		//se cerco di settare il nuovo start su una mask continua a provare

		tmp.setMotifWidth(); //setta il motifWidth guardando alla prima sequenza inserita
		//rimuovi la sequenza presa dallo stack di input se tutto e' andato bene
		inputSequences.remove(startingMotif);

		int shortestLength = 9999999;

		if(startingNo < 1){
			startingNo = 1;
		}
		for (int i = 0; i < (int)startingNo ; i++) {//da decidere
			if (Debug.ON) System.out.println("adding sequence " + i);
			Motif seq = new Motif();
			seq = inputSequences.get(0); //get Motif i
			maskGoFlag=false;
			int start = 0;

			while(!maskGoFlag){
				start=BEARManager.searchInitialMotif(startingMotif.extractMotifFromSequence(), seq, m, motifLength, true);
				if (start == -1){
					break;
				}
				for(int j=start; j < start+motifLength; j++){
					maskGoFlag= !seq.isMasked(j); //mGF rimane diventa true per ogni posizione corretta. Appena diventa false, viene catchata e riparte il while
					if (Debug.ON) System.out.println(start + "\t" + motifLength + "\t" + maskGoFlag + "\t"+ seq.getMotifEnd()+ "\t" + seq.getMotifStart()+ "\t" + seq.printMask());
					if (!maskGoFlag) break;
				}

			}
			if(start+motifLength > seq.getSequence().length()) System.err.println("tentato motifStart maggiore del massimo consentito: " 
					+ start + ", " + "lunghezza seq: " + seq.getSequence().length());


			if (start != -1){ //fai questa serie di operazioni solo se si puo', altrimenti ignora la sequenza.
				seq.setMotifStart(start);
				seq.setMotifEnd(start+motifLength, motifLength);
				if( !Debug.checkStartEnd_element(seq, motifLength)){
					System.err.println(seq.getName() + "\t" + seq.getMotifStart() + "\t"+ seq.getMotifEnd() + "\t"+ start + "\t"+ motifLength + "\t"+ tmp.getMotifWidth());
					seq.setMotifEnd(start+motifLength);

				}
				tmp.addMotif(seq);
				int tmpShortest = seq.getSequence().length();
				if (tmpShortest < shortestLength){
					shortestLength = tmpShortest;
				}
				if (tmp.getListMotif().get(tmp.cardinality()-1).getMotifEnd() - tmp.getListMotif().get(tmp.cardinality()-1).getMotifStart() != motifLength ){
					System.err.println("errore inizializzazione");
					tmp.removeMotif(tmp.cardinality()-1);
				}else{
					inputSequences.remove(seq);

				}
			}
		}
		tmp.setShortestSequenceLength(shortestLength);
		return tmp;
	}

	public static int getRandomIndex(MotifHandler m){
		return (int)(Math.random()*m.cardinality());
	}


	public static int getPoisson(double lambda) {
		double L = Math.exp(-lambda);
		double p = 1.0;
		int k = 0;
		do {
			k++;
			p *= Math.random();
		} while (p > L);
		return k - 1;
	}


	private static boolean addSequenceToMotifWithPSSM(MotifHandler mh, ArrayList<Motif> inputSequences, double[][] mat){
		Motif m = new Motif();
		//FIX__7/02__ per i dataset disomogenei, per evitare di prendere seq pi� corta del motivo
		int count = 0;
		int triesLim = 100;
		do{
			Collections.shuffle(inputSequences, new Random(System.currentTimeMillis()));
			//System.out.println("inputSequences size: " + inputSequences.size());
			m = inputSequences.get(0);
			count++;
		}while( (m.getSequence().length() < mh.getMotifWidth() ) && count <= triesLim);

		//System.out.println(m.getName());

		if(count > triesLim){
			if(mh.cardinality() == 0){
				System.err.println("Non sono riuscito a trovare una sequenza lunga almeno " + mh.getMotifWidth() + ". Inoltre tutte quelle di prima facevano cagare quindi esco. Il motivo non c'�. (Doppio senso, pero' a questo punto il motivo c'� se sono uscito. PARADOX)");
				System.exit(-1);
			}else{
				return false;
			}
			//			System.out.println("DELTA2 = " + PSSM.computeDeltaScoreAdd(mh.getMotifWidth(), PSSM.computeCounts(mh), alphas, dataPriors));

		}

		//seq.clearStartEndIndexes();

		PSSM.computePSSM(mh, inputSequences);

		boolean maskGoFlag=false;
		int res = 0;
		int motifLength = mh.getMotifWidth();
		while(!maskGoFlag){
			//estrai posizione ottimale da pssm, la funzione da punteggio -inf alle pos mascherate
			res=BEARManager.searchMotifUsingPSSM(mh.getPSSM(), m, mat ,motifLength, true);
			if (res==-1){
				break;
			}
			for(int j=res; j <= res+motifLength; j++){ //NON DEVE ESSERE NECESSARIO
				maskGoFlag= !m.isMasked(j);//mGF diventa true per ogni posizione corretta. Appena diventa false, viene catchata e riparte il while
				//System.out.println(maskGoFlag + "\t"+ m.getMotifEnd()+ "\t" + m.getMotifStart()+ "\t" + m.printMask());
				if (!maskGoFlag) break;
			}
		}
		if(res!=-1){
			m.setMotifStart(res, motifLength);
			m.setMotifEnd(res+motifLength, motifLength);

			//System.out.println("sequenza aggiunta " + seq.toString() + " " + res);

			//Se la nuova sequenza e' piu' corta dell'attuale piu' corta, allora aggiorna shortestSequenceLength, SI, A MANO.
			int tmp = m.getSequence().length();
			if (tmp < mh.getShortestSequenceLength()){
				mh.setShortestSequenceLengthPrev(mh.getShortestSequenceLength());
				mh.setShortestSequenceLength(tmp);
			}		

			inputSequences.remove(m);
			mh.addMotif(m);
			return true;
		}else{
			return false;
		}
	}

	private static void undoAdd(MotifHandler mh, ArrayList<Motif> inputSequences){
		inputSequences.add(mh.getObjectMotif(mh.cardinality()-1));
		mh.removeMotif(mh.cardinality()-1);

	}

	private static void removeSequenceFromMotifTestPSSM(MotifHandler mh, ArrayList<Motif> inputSequences){
		int rand=getRandomIndex(mh);
		inputSequences.add(mh.getObjectMotif(rand));
		//System.out.println(mh.getObjectMotif(rand));
		mh.removeMotif(rand);
	}

	//	private static void removeSequenceFromMotifTestPSSM(MotifHandler mh, ArrayList<Motif> inputSequences){
	//TEST PER RIMUOVERE SOLO LE ULTIME
	//		inputSequences.add(mh.getObjectMotif(mh.cardinality()-1));
	//		System.out.println(mh.getObjectMotif(mh.cardinality()-1));
	//		mh.removeMotif(mh.cardinality()-1);
	//	}

	private static void undoRemove(MotifHandler mh, ArrayList<Motif> inputSequences){
		mh.addMotif(inputSequences.get(inputSequences.size()-1));
		inputSequences.remove(inputSequences.size()-1);
	}



	private static void undoChangedIndex(MotifHandler mh) {
		for(Motif m:mh.getListMotif()){
			m.setMotifStartUndo();
			m.setMotifEndUndo();
		}
		mh.setMotifWidth(mh.getMotifWidthPrev());

	}

	private static void changeSubopt(MotifHandler mh, ArrayList<Motif> inputSequences, double[][] mat ){


		//Scegli sequenza a caso
		int rand=getRandomIndex(mh);
		//spostala indietro nelle input
		inputSequences.add(mh.getObjectMotif(rand));
		//toglila da mh
		mh.removeMotif(rand);
		//risetta che non se sa mai
		mh.setMotifWidth();

		//cosi' la pssm non e' influenzata dalla subopt tolta //OTTIMIZZABILE
		PSSM.computePSSM(mh, inputSequences); 

		mh.addMotif(inputSequences.get(inputSequences.size()-1)); //riprendo la sequenza
		inputSequences.remove(inputSequences.size()-1); //e la rimuovo dalle input



		boolean ok = false;
		int count = 0;

		//1 tentativi o condizione di uscita
		while(!ok && count < 1){

			int newIndex = (int)(Math.random() * mh.getObjectMotif(mh.cardinality()-1).getSequenceList().size());

			//se la sequenza e' abbastanza lunga accetta senn� torna indietro
			if(mh.getObjectMotif(mh.cardinality()-1).getSequenceList().get(newIndex).length() > mh.getMotifWidth() ){

				//Setta il nuovo indice alla sequenza ripresa
				mh.getObjectMotif(mh.cardinality()-1).setIndex(newIndex);

				ok=true;

				boolean maskGoFlag=false;
				int res = 0;
				int motifLength = mh.getMotifWidth();

				while(!maskGoFlag){

					res=BEARManager.searchMotifUsingPSSM(mh.getPSSM(), mh.getObjectMotif(mh.cardinality()-1), mat, mh.getMotifWidth());

					for(int j=res; j <= res+motifLength; j++){
						maskGoFlag= !mh.getObjectMotif(mh.cardinality()-1).isMasked(j); //mGF rimane diventa true per ogni posizione corretta. Appena diventa false, viene catchata e riparte il while
						if (!maskGoFlag) break;
					}
				}



				res=BEARManager.searchMotifUsingPSSM(mh.getPSSM(), mh.getObjectMotif(mh.cardinality()-1), mat, motifLength);
				if (res != -1){
					mh.getObjectMotif(mh.cardinality()-1).setMotifStart(res, motifLength);
					mh.getObjectMotif(mh.cardinality()-1).setMotifEnd(res+motifLength, motifLength);
					ok = false;
					count++;
				}
				//System.out.println();

			}else{
				count++;
			}
		}

		if (ok == false){
			mh.getObjectMotif(mh.cardinality()-1).setIndexUndo();
			mh.getObjectMotif(mh.cardinality()-1).setMotifStart(mh.getObjectMotif(mh.cardinality()-1).getMotifStart());
			mh.getObjectMotif(mh.cardinality()-1).setMotifEnd(mh.getObjectMotif(mh.cardinality()-1).getMotifEnd());
			return;
		}

		//ATTENZIONE A QUANDO LE SUBOPT SONO DIVERSE STRUTTURE LOCALI DI GRANDEZZA DIVERSA
		//aggiungendo le maschere non � pi� possibile questa situazione


	}

	private static void undoChangeSubopt(MotifHandler mh) {
		if (Debug.VERBOSE) System.out.println("subopt UNDO -------");

		mh.getObjectMotif(mh.cardinality()-1).setIndexUndo();

		mh.getObjectMotif(mh.cardinality()-1).setMotifStartUndo();
		mh.getObjectMotif(mh.cardinality()-1).setMotifEndUndo();


	}

	private static void undoPartials(MotifHandler mh){
		if (Debug.VERBOSE) System.out.println("partial UNDO ------");
		for (Motif m: mh.getListMotif()){
			m.setPartialUndo();
		}
		return;
	}

	public static void ctrlZ(int operazione,MotifHandler mh, ArrayList<Motif> inputSequences){
		if(operazione==0){
			undoAdd(mh,inputSequences);
		}else if(operazione==1){
			undoRemove(mh,inputSequences);
		}else if(operazione==6){
			undoChangeSubopt(mh);
		}else{ //casi 2,3,4,5 e' solo un cambio indici (stretch <- -> -><- shift e recalculate)
			undoChangedIndex(mh);
		}
		undoPartials(mh);

	}

	private static void recalculateMotif(MotifHandler mh, double[][] mat){
		int choosenOne = getRandomIndex(mh);
		int motifWidth = mh.getMotifWidth();
		String tmpMotif=mh.getSequenceMotif(choosenOne);
		String tmpName=mh.getObjectMotif(choosenOne).getName();
		mh.getObjectMotif(choosenOne).setMotifStart(mh.getObjectMotif(choosenOne).getMotifStart());//questi due servono per non far sballare la macchina
		mh.getObjectMotif(choosenOne).setMotifEnd(mh.getObjectMotif(choosenOne).getMotifEnd());

		//System.out.println("lunghezza motivo seq " + (mh.getListMotif().get(choosenOne).getMotifEnd() - mh.getListMotif().get(choosenOne).getMotifStart()));
		//System.out.println("lugnhezza motivo handler " + mh.getMotifLength());
		for (Motif element: mh.getListMotif()){
			//System.out.println(element.getMotif());
			if (!(element.getName().equals(tmpName))){

				int newMotifStart = BEARManager.searchInitialMotif(tmpMotif, element, mat,  motifWidth);
				//System.out.println(newMotifStart);
				element.setMotifStart(newMotifStart, motifWidth);
				element.setMotifEnd(newMotifStart + motifWidth, motifWidth);
				//				if(!element.setMotifEnd(newMotifStart + motifWidth)){
				//					System.err.println((newMotifStart + motifWidth) + "\t" + element.printMask());}
			}
		}
	}

	private static void shiftWindow(MotifHandler mh){
		//shifta tutte le finestre a sinistra o a destra
		int shiftDirection = +1;
		int motifWidth=mh.getMotifWidth();
		if (Math.random()<0.5){
			shiftDirection = -1;
		}
		int shiftIntensity = getPoisson(1);
		for(Motif m: mh.getListMotif()){
			int sh = shiftIntensity;
			if (m.getMotifStart()+sh*shiftDirection < 0 || m.getMotifEnd()+sh*shiftDirection > m.getSequence().length()){
				sh=0; //se esce dai bordi
			}else{
				for (int j=1; j<=shiftIntensity; j++){
					if (m.isMasked(m.getMotifStart()+j*shiftDirection) || m.isMasked(m.getMotifEnd()+j*shiftDirection)){
						sh = 0;
						//se una finestra non puo' essere spostata, rimane ferma
					}
				}
			}

			shiftSingleMotif(m,sh,shiftDirection, motifWidth);
			if (!Debug.checkStartEnd_element(m, mh.getMotifWidth()) ){
				System.err.println(m.getName() + "\t"
						+ "shiftato male!");
			}
		}

	}


	private static void shiftSingleMotif(Motif m, int shift, int shiftDirection, int motifWidth){
		m.setMotifStart(m.getMotifStart()+shift*shiftDirection, motifWidth);
		m.setMotifEnd(m.getMotifEnd()+shift*shiftDirection, motifWidth);
	}

	private static void expand(MotifHandler mh, int widthUpperBound){
		int resize = 1;
		resize = getPoisson(1); //e' la media della distribuzione di poisson
		boolean fromTheRight = true;
		int motifWidth=mh.getMotifWidth();

		if (Math.random() < 0.5){
			fromTheRight = false;
		}
		if( mh.getMotifWidth() + resize > widthUpperBound){
			resize = 0;
		}
		//		if (expand){
		if (Debug.ON) System.out.println("try expand " + resize);
		boolean doItFlag = true;
		for (Motif element: mh.getListMotif()){ //controlla se e' possibile espandere tutte le sequenze altrimenti nada
			if (fromTheRight && element.getMotifEnd()+resize > element.getSequence().length()){
				doItFlag = false;
			}else if(!fromTheRight && element.getMotifStart() - resize < 0){
				doItFlag = false;
			}else if(fromTheRight && element.isMasked(element.getMotifEnd() + resize)){
				doItFlag = false;
			}else if(!fromTheRight && element.isMasked(element.getMotifStart() - resize)){
				doItFlag = false;
			}
		}
		//			System.out.println(doItFlag);
		if(!doItFlag){
			resize=0;
			if (Debug.ON) System.out.println("but no expand");

		}
		for (Motif element: mh.getListMotif()){

			//Probabilmente qua c'� un problema se prova ad andare su una maschera, perch� in caso alcuni singoli non li fa.
			if ((fromTheRight && !((element.getMotifEnd() + resize) > (element.getSequence().length()))) || (!fromTheRight && ((element.getMotifStart() - resize) < 0)  )){
				element.setMotifEnd(element.getMotifEnd() + resize , motifWidth);
				element.setMotifStart(element.getMotifStart(), motifWidth);
			}else if((!fromTheRight && !((element.getMotifStart() - resize) < 0)) || (fromTheRight && ((element.getMotifEnd() + resize) > element.getSequence().length()))){
				element.setMotifStart(element.getMotifStart() - resize, motifWidth);
				element.setMotifEnd(element.getMotifEnd(), motifWidth);
			}
		}
		mh.setMotifWidth();
	}



	private static void shrink(MotifHandler mh, int widthLowerBound){
		//SHRINK NON OMOGENEO: mh prendeva piu' volte le stesse sequenze da inputsequences, problema risolto... se magari
		int resize = 1;
		resize = getPoisson(1); //2 e' la media della distribuzione di poisson
		boolean fromTheRight = true;
		int motifWidth=mh.getMotifWidth();
		if (Math.random() < 0.5){
			fromTheRight = false;
		}

		if( mh.getMotifWidth() - resize < widthLowerBound){
			resize = 0;
		}

		if (Debug.ON) System.out.println("shrink: " + resize);
		boolean ctrlz=false;
		for (Motif element: mh.getListMotif()){	


			if (fromTheRight){
				if (!element.setMotifEnd(element.getMotifEnd() - resize, motifWidth)){
					ctrlz=true;

				}
				element.setMotifStart(element.getMotifStart(), motifWidth);
			}else{
				if (!element.setMotifStart(element.getMotifStart() + resize, motifWidth)){
					ctrlz=true;

				}
				element.setMotifEnd(element.getMotifEnd(), motifWidth);
			}

			//				System.out.println("Dopo: " + element.extractMotifFromSequence());


			//		System.out.println(mh.getMotifWidth());

		}
		if (ctrlz){

			System.err.println("reverting indexes after shrink");
			undoChangedIndex(mh);
		}
		mh.setMotifWidth();
		//		System.out.println(mh.getMotifWidth());
	}


	public int perturbateMotif(MotifHandler mh, ArrayList<Motif> inputSequences, double[][] mat, int widthUpperBound, boolean weHaveSubopt, int model_limit ) {
		//questo chiama le funzioni di sopra
		/*
		 * 0:Add
		 * 1:Remove
		 * 2:Shift,Enlarge,Shrink
		 */
		int counta = 0;
		int widthLowerBound = 3;
		int minSeq = 2;
		int ops = 6;
		if (weHaveSubopt) ops= 7;
		Random generator = new Random();

		if (mh.cardinality() > minSeq && mh.cardinality() < model_limit && inputSequences.size() >0){
			int roll = generator.nextInt(ops);
			switch(roll){
			case 0: 
				if (Debug.VERBOSE) System.out.println("---Sequence Removed---");
				removeSequenceFromMotifTestPSSM(mh,inputSequences);
				return 1;
			case 1: 
				if (Debug.VERBOSE) System.out.println("---Sequence Adding---");
				addSequenceToMotifWithPSSM(mh, inputSequences, mat);
				if (Debug.VERBOSE) System.out.println("---Sequence Added---");
				return 0;
			case 2:	
				if (Debug.VERBOSE) System.out.println("---Shifting Window---");
				shiftWindow(mh);
				return 2;
			case 3: 
				if (Debug.VERBOSE) System.out.println("---Recalculating Motif---");
				recalculateMotif(mh,mat);
				return 3;
			case 4: 
				if (Debug.VERBOSE) System.out.println("---Expanding Motif---");
				expand(mh, widthUpperBound);
				return 4;
			case 5: 
				if (Debug.VERBOSE) System.out.println("---Shrinking Motif---");
				shrink(mh,widthLowerBound);
				return 5;
			case 6: 
				if (Debug.VERBOSE) System.out.println("---Changing Subopt---");
				changeSubopt(mh, inputSequences, mat);
				return 6;
				//				case 6: System.out.println("---Jumping---");
				//				jumpOnNewMotif(mh,mat);
				//				return 5;
			default: return -1;
			}
		}else if (mh.cardinality() <= minSeq){
			if (Debug.VERBOSE) System.out.println("---Forcing Sequence Added---");
			while(counta < 100 && !addSequenceToMotifWithPSSM(mh, inputSequences, mat)){
				counta++;
			}
			if (counta==100){
				return -1;
			}
			if (Debug.VERBOSE) System.out.println("---Forcing Seq Added Uscita---");
			//System.out.println(counta);

			return 0;
		}else{
			int roll = generator.nextInt(ops-1);
			switch(roll){
			case 0: 
				if (Debug.VERBOSE) System.out.println("---Sequence Removed---");
				removeSequenceFromMotifTestPSSM(mh,inputSequences);
				return 1;
			case 1:	
				if (Debug.VERBOSE) System.out.println("---Shifting Window---");
				shiftWindow(mh);
				return 2;
			case 2: 
				if (Debug.VERBOSE) System.out.println("---Recalculating Motif---");
				recalculateMotif(mh,mat);
				return 3;
			case 3: 
				if (Debug.VERBOSE) System.out.println("---Expanding Motif---");
				expand(mh,widthUpperBound);
				return 4;
			case 4: 
				if (Debug.VERBOSE) System.out.println("---Shrinking Motif---");
				shrink(mh,widthLowerBound);
				return 5;
			case 5: 
				if (Debug.VERBOSE) System.out.println("---Changing Subopt---");
				changeSubopt(mh, inputSequences, mat);
				return 6;
				//				case 5: System.out.println("---Jumping---");
				//				jumpOnNewMotif(mh,mat);
				//				return 5;
			default: return -1;
			}
		}
	}



	/*











	 */
}
