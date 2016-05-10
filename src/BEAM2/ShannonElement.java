package BEAM2;

public class ShannonElement {
	char name;
	double prior;
	
	ShannonElement(char name_, double prior_){
		this.name = name_;
		this.prior = prior_;
	}

	public void setPrior(double d) {
		this.prior = d;
	}
	public double getPrior(){
		return prior;
	}

	public char getName() {
		return name;
	}

	public void setName(char name) {
		this.name = name;
	}
}
