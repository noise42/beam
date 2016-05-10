package BEAM2;

import java.util.ArrayList;

public class Shannon {
	private ArrayList<ShannonElement> priorsVector;
//	private double[] pseudo = new double[83];
	
	Shannon(){
		String bearCharacters = "abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";
		this.priorsVector = new ArrayList<ShannonElement>();
		for(int i=0; i<83;i++){
			priorsVector.add(new ShannonElement(bearCharacters.charAt(i), 1.0/83 ) );
		}
	}
	
	public void fillFromInput(ArrayList<Motif> inputSequences){
		String bearCharacters = "abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";

		int[] count = new int[83];
		for(Motif m: inputSequences){
			char[] chars = m.getSequence().toCharArray();
			for(char c : chars){
				count[bearCharacters.indexOf(c)] += 1;

			}
		}
		int total=0;
		for(int i=0; i<83;i++){
			total += count[i];
		}
		for(int i=0; i<83;i++){
			priorsVector.get(i).setPrior((double)count[i]/total);
		}
	}
	
	public void fillFromFile(double [] dataPriors){
		for(int i=0; i<83;i++){
			priorsVector.get(i).setPrior(dataPriors[i]);
		}
	}
	
	
	@Override
	public String toString(){
		String tmp = "";
		double total=0;
		for(int i = 0 ; i<priorsVector.size(); i++){
			tmp+= priorsVector.get(i).getName() + "\t" +priorsVector.get(i).getPrior() + "\n";
			total += priorsVector.get(i).getPrior();
		}
		tmp+="\nTotale " + total;
		return tmp;
	}

	public ArrayList<ShannonElement> getPriorsVector() {
		return priorsVector;
	}
	
	public double getPriorFromChar(char c){
		String bearCharacters = "abcdefghi=lmnopqrstuvwxyz^!\"#$%&\'()+234567890>[]:ABCDEFGHIJKLMNOPQRSTUVW{YZ~?_|/\\}@";

		return priorsVector.get(bearCharacters.indexOf(c)).getPrior();
	}

	public void setPriorsVector(ArrayList<ShannonElement> priorsVector) {
		this.priorsVector = priorsVector;
	}
}
