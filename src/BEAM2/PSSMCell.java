package BEAM2;

//import java.util.ArrayList;
/*un istanza PSSM e' una cella della PSSM*/
public class PSSMCell {
	char name;
	double occurrence=0.0;
	
	public PSSMCell(char name, double i) {
		this.name=name;
		this.occurrence=i;
	}
}
