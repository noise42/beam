package BEAM2;


public class MotifMemory {
	private MotifHandler handlerMemory;
	private double score;
	
	public MotifMemory(){
		this.handlerMemory = new MotifHandler();
		this.score = 0.0;
	}
	
	public MotifMemory(double score_, MotifHandler mh_) {
		this.handlerMemory = new MotifHandler();
		this.handlerMemory.MotifHandlerClone(mh_);
		
		this.score = score_;
	}

	public void finalize(){
		
	}
	
	public boolean tryMask(double score_, MotifHandler mh_){
		//se lo score � maggiore di quello in memoria, salva quell'mh
		if(score_ > this.score){
			this.handlerMemory.MotifHandlerClone(mh_);
			this.score = score_;
			return true;
		}else{
			return false;
		}
	}
/*	
	public void applyMask(ArrayList<Motif> inputSequences, String name, int start, int end){
		
	}
*/
	public MotifHandler getHandlerMemory() {
		return handlerMemory;
	}

	public void setHandlerMemory(MotifHandler handlerMemory) {
		this.handlerMemory = handlerMemory;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
}
