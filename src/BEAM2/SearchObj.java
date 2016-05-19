package BEAM2;
import java.util.ArrayList;
import java.util.*;

//Contiene la struttura dati necessaria al search after run, con pssm e dati di confusion matrix
public class SearchObj {
	private ArrayList<ArrayList<PSSMCell>> PSSMdata;
	private ArrayList<ArrayList<PSSMCell>> qpssm;

	private int TP,FP,TN,FN;
	private double mediaPart;
	private double pvalue;
	
	SearchObj(ArrayList<ArrayList<PSSMCell>> pssmdata, int tp, int fp, int tn, int fn, double mediapart){
		//copio uno ad uno i valori di pssmdata (altrimenti si porta appresso la referenza
		//java di merda 2.0
		this.qpssm = new ArrayList<ArrayList<PSSMCell>>();
		for (ArrayList<PSSMCell> qaL : qpssm) {
			ArrayList<PSSMCell> tmpqAL = new ArrayList<PSSMCell>();
			for (PSSMCell cell : qaL){
				tmpqAL.add(cell);
			}
		}
		this.PSSMdata = new ArrayList<ArrayList<PSSMCell>>();
		for (ArrayList<PSSMCell> aL : pssmdata) {
			ArrayList<PSSMCell> tmpAL = new ArrayList<PSSMCell>();
			for (PSSMCell cell : aL){
				tmpAL.add(cell);
			}
			this.PSSMdata.add(tmpAL);
			//TODO non sono affatto sicuro che questo workaround funzioni
		}
		
		//copio i valori di confusion matrix (easy)
		this.TP = tp;
		this.FP = fp;
		this.TN = tn;
		this.FN = fn;
		
		this.mediaPart=mediapart;
		
		this.pvalue=-1;
	}
	
	public double computeMCC(){
		//calcola MCC avendo i propri valori di tp fp tn fn
		//TP non pu� essere 0 perch� sono al minimo i costituenti dell'MH
		//FP pu� essere 0 in casi molto positivi
		//TN pu� essere 0 in casi molto negativi
		//FN pu� essere 0 la maggior parte delle volte
		//pseudoconte:
		if (FP==0) FP++;
		if (TN==0) TN++;
		if (FN==0) FN++;
		//return (TP*TN - FP*FN)/Math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
		return (Math.exp(Math.log(TP)+Math.log(TN)) - Math.exp(Math.log(FP)+Math.log(FN)))/Math.exp(.5*(Math.log(TP+FP)+Math.log(TP+FN)+Math.log(TN+FP)+Math.log(TN+FN)));
	}
	
	public String computeConsensus() {
		//legge PSSMdata e prende la consensus, se ci sono frazioni SIMILI e ALTE
		//scrive entrambi. 
		//SIMILI: stessa %?
		//ALTE: oltre 35%?
		//TODO ^^^
		char c, extraC;
		double freq;
		
		String consensus="";
		//int column = 0;
		//creo una lista di indici in cui segnare i char della seconda riga
	//	ArrayList<Integer> extraIdx = new ArrayList<Integer>();
		String secondRow="";
		StringBuilder builder = new StringBuilder(secondRow);

		for (ArrayList<PSSMCell> AL: PSSMdata){
			freq=0.0;
			c = ' ';
			extraC = ' ';
			//per ogni colonna dell'allineamento determina quella a % + alta
			for (PSSMCell cell: AL){
				//trovo il massimo
				if (cell.occurrence >= freq){
					freq = cell.occurrence;
					//gestire caso in cui trova prima un 41% e poi un 45% (ad es)
					//metto il precedente nel secchio e via
					extraC = c; 
					c = cell.name;
				}else if (cell.occurrence >=0.35){
					//se � maggiore del 40% ma minore del massimo trovato comunque segna e metti su seconda riga
					extraC = cell.name;
				}
			}
			if (freq < 0.35){
				consensus += '*';
			}else{
				consensus += c;
			}
			
			builder.append(extraC);
			
		}
		return consensus + "\n" + secondRow + "\n";
	}
	
	public String computeqbearConsensus() {
		//legge PSSMdata e prende la consensus, se ci sono frazioni SIMILI e ALTE
		//scrive entrambi. 
		//SIMILI: stessa %?
		//ALTE: oltre 35%?
		//TODO ^^^
		char c, extraC;
		double freq;
		ArrayList<ArrayList<PSSMCell>> tmpqpssm = new ArrayList<ArrayList<PSSMCell>>();
		for (ArrayList<PSSMCell> AL: PSSMdata){
			ArrayList<PSSMCell> qrow = new ArrayList<PSSMCell>();
			for (PSSMCell pcell: AL){
				qrow.add(new PSSMCell(qBEAR.identifyClass(pcell.name), pcell.occurrence));
			}
			tmpqpssm.add(qrow);
			

		}
		Map<Character,Double> dict = new HashMap<Character,Double>();
		for(ArrayList<PSSMCell> AL: tmpqpssm){
			dict.clear();
			for(PSSMCell pcell: AL){
				
				if(dict.get(pcell.name)==null){
					dict.put(pcell.name, pcell.occurrence);
				}else{
					//double value = dict.get(pcell.name);
					dict.put(pcell.name, dict.get(pcell.name) + pcell.occurrence);
				}
			}
			ArrayList<PSSMCell> qrow = new ArrayList<PSSMCell>();
			for (Character k: dict.keySet()){
				qrow.add(new PSSMCell(k, dict.get(k)) );
			}
			qpssm.add(qrow);
			
			/*
			for(PSSMCell pc : qrow){
				System.out.println(pc.name + ":" + pc.occurrence);
			}
			System.out.println("fine riga");
			*/
			
		}
		
		String consensus="";
		//creo una lista di indici in cui segnare i char della seconda riga
		String secondRow="";
		StringBuilder builder = new StringBuilder(secondRow);

		for (ArrayList<PSSMCell> AL: qpssm){
		
			freq=0.0;
			c = ' ';
			extraC = ' ';
			//per ogni colonna dell'allineamento determina quella a % + alta
			for (PSSMCell cell: AL){
				//trovo il massimo
				if (cell.occurrence >= freq){
					freq = cell.occurrence;
					//gestire caso in cui trova prima un 41% e poi un 45% (ad es)
					//metto il precedente nel secchio e via
					extraC = c; 
					c = cell.name;
				}else if (cell.occurrence >=0.35){
					//se � maggiore del 40% ma minore del massimo trovato comunque segna e metti su seconda riga
					extraC = cell.name;
				}
			}
			if (freq < 0.35){
				consensus += '*';
			}else{
				consensus += c;
			}
			
			builder.append(extraC);
			
		}
		return consensus + "\n" + secondRow + "\n";
		
	}
	
	int getTP(){
		return this.TP;
	}
	int getFP(){
		return this.FP;
	}
	int getTN(){
		return this.TN;
	}
	int getFN(){
		return this.FN;
	}
	
	double getMeanPart(){
		return this.mediaPart;
	}

	public void setPvalue(double pvalue_) {
		this.pvalue = pvalue_;
	}
	
	public double getPvalue(){
		return this.pvalue;
	}
	




}
