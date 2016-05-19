package BEAM2;

import java.util.ArrayList;

public class Debug {
	private Debug() {}

	public static final boolean ON = false;
	//	public static final boolean ON = true;


	public static boolean VERBOSE = false;
	
	public static void check(int idx){
		System.out.println("check " + idx);
	}

	static public void setVERBOSE(boolean state){
		VERBOSE = state;
	}

	static public ArrayList<Motif> checkStartEnd(MotifHandler mh, ArrayList<Motif> inputSequences){
		ArrayList<Motif> list = new ArrayList<Motif>();

		for(Motif m: mh.getListMotif()){
			if(m.getMotifEnd()-m.getMotifStart() != mh.getMotifWidth()){
				list.add(m);
			}
		}
		//sposta le incriminate nelle input
		for (Motif move: list){
			inputSequences.add(move);
			mh.getListMotif().remove(move);
		}
		return list;
	}
	static public boolean checkPreviousStartEnd(MotifHandler mh){

		for(Motif m: mh.getListMotif()){
			if(m.getMotifEndPrev()-m.getMotifStartPrev() != mh.getMotifWidthPrev()){
				return false;
			}
		}

		return true;
	}

	static public ArrayList<Motif> checkMasks(MotifHandler mh, ArrayList<Motif> inputSequences){
		ArrayList<Motif> list = new ArrayList<Motif>();
		for(Motif m: mh.getListMotif()){
			for (int j=m.getMotifStart(); j<=m.getMotifEnd(); j++){
				if (m.isMasked(j)){
					list.add(m);

				}
			}

		}
		for (Motif move: list){
			inputSequences.add(move);
			mh.getListMotif().remove(move);
		}
		return list;
	}
	public static boolean checkMask_element(Motif element) {
		for (int j=element.getMotifStart(); j<=element.getMotifEnd(); j++){
			if (element.isMasked(j)){
				return false;
			}
			return false;
		}
		return true;
	}
	public static boolean checkStartEnd_element(Motif m, int motifWidth){
		return (m.getMotifEnd()-m.getMotifStart()) == motifWidth; 
	}

	public static void checkRNA(MotifHandler mh) {
		// controlla se ci sono due rna con lo stesso nome
		ArrayList<String> rna = new ArrayList<String>();
		for (Motif m: mh.getListMotif()){
			if (rna.contains( m.getName() )  ){
				System.out.println(mh.toString());
				IO.scream(-1);
			}else{
				rna.add(m.getName());
			}
		}

	}
}
