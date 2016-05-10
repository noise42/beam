package BEAM2;

public class qBEAR {
	static public char identifyClass(char BEARchar){
//		int index = 0;
		String stemClassShort = "abcde";
		String stemClassMedium = "fghi";
		String stemClassLong = "=";

		String loopClassShort = "jklmnopqr";
		String loopClassMedium = "stuvwxyz";
		String loopClassLong = "^";

		String internalLoopClassShort = "!\"#$%23456";
		String internalLoopClassMedium = "&\'()7890";
		String internalLoopClassLong = "+>";

		String bulgeClass = "[]";
		
		String stemBranchClassShort = "ABCDE";
		String stemBranchClassMedium = "FGHI";
		String stemBranchClassLong = "J";

		String internalLoopBranchClassShort = "KLMNYZ~?";
		String internalLoopBranchClassMedium = "OPQRS_|/\\";
		String internalLoopBranchClassLong = "TUVWYZ@";

		String bulgeBranchClass = "{}";
		
		String branchClass = ":";
		
		if (stemClassShort.contains(Character.toString(BEARchar)) ){
			return 'z';
		}else if(stemClassMedium.contains(Character.toString(BEARchar)) ){
			return 'a';
		}else if(stemClassLong.contains(Character.toString(BEARchar)) ){
			return 'q';
			
		}else if(loopClassShort.contains(Character.toString(BEARchar)) ){
			return 'x';
		}else if(loopClassMedium.contains(Character.toString(BEARchar)) ){
			return 's';
		}else if(loopClassLong.contains(Character.toString(BEARchar)) ){
			return 'w';
			
		}else if(internalLoopClassShort.contains(Character.toString(BEARchar)) ){
			return 'c';
		}else if(internalLoopClassMedium.contains(Character.toString(BEARchar)) ){
			return 'd';
		}else if(internalLoopClassLong.contains(Character.toString(BEARchar)) ){
			return 'e';
			
		}else if(bulgeClass.contains(Character.toString(BEARchar)) ){
			return 'b';
			
		}else if(stemBranchClassShort.contains(Character.toString(BEARchar)) ){
			return 'v';
		}else if(stemBranchClassMedium.contains(Character.toString(BEARchar)) ){
			return 'f';
		}else if(stemBranchClassLong.contains(Character.toString(BEARchar)) ){
			return 'r';
			
		}else if(internalLoopBranchClassShort.contains(Character.toString(BEARchar)) ){
			return 'n';
		}else if(internalLoopBranchClassMedium.contains(Character.toString(BEARchar)) ){
			return 'h';
		}else if(internalLoopBranchClassLong.contains(Character.toString(BEARchar)) ){
			return 'y';
			
		}else if(bulgeBranchClass.contains(Character.toString(BEARchar)) ){
			return 'g';
			
		}else if(branchClass.contains(Character.toString(BEARchar)) ){
			return 't';
		}else{
			return '-';
		}
	}
	
	static public char partialIdentifyClass(char BEARchar){
//		int index = 0;
		String stemClassShort = "abcde";
		String stemClassMedium = "fghi";
		String stemClassLong = "=";

		String loopClassShort = "jklmnopqr";
		String loopClassMedium = "stuvwxyz";
		String loopClassLong = "^";

		String internalLoopClassShort = "!\"#$%23456";
		String internalLoopClassMedium = "&\'()7890";
		String internalLoopClassLong = "+>";

		String bulgeClass = "[]";
		
		String stemBranchClassShort = "ABCDE";
		String stemBranchClassMedium = "FGHI";
		String stemBranchClassLong = "J";

		String internalLoopBranchClassShort = "KLMNYZ~?";
		String internalLoopBranchClassMedium = "OPQRS_|/\\";
		String internalLoopBranchClassLong = "TUVWYZ@";

		String bulgeBranchClass = "{}";
		
		String branchClass = ":";
		
		if (stemClassShort.contains(Character.toString(BEARchar)) ){
			return BEARchar;
		}else if(stemClassMedium.contains(Character.toString(BEARchar)) ){
			return BEARchar;
		}else if(stemClassLong.contains(Character.toString(BEARchar)) ){
			return BEARchar;
			
		}else if(loopClassShort.contains(Character.toString(BEARchar)) ){
			return 'x';
		}else if(loopClassMedium.contains(Character.toString(BEARchar)) ){
			return 's';
		}else if(loopClassLong.contains(Character.toString(BEARchar)) ){
			return 'w';
			//TODO Finire pBEAR, senza floyd non ricordo 
		}else if(internalLoopClassShort.contains(Character.toString(BEARchar)) ){
			return 'c';
		}else if(internalLoopClassMedium.contains(Character.toString(BEARchar)) ){
			return 'd';
		}else if(internalLoopClassLong.contains(Character.toString(BEARchar)) ){
			return 'e';
			
		}else if(bulgeClass.contains(Character.toString(BEARchar)) ){
			return 'b';
			
		}else if(stemBranchClassShort.contains(Character.toString(BEARchar)) ){
			return 'v';
		}else if(stemBranchClassMedium.contains(Character.toString(BEARchar)) ){
			return 'f';
		}else if(stemBranchClassLong.contains(Character.toString(BEARchar)) ){
			return 'r';
			
		}else if(internalLoopBranchClassShort.contains(Character.toString(BEARchar)) ){
			return 'n';
		}else if(internalLoopBranchClassMedium.contains(Character.toString(BEARchar)) ){
			return 'h';
		}else if(internalLoopBranchClassLong.contains(Character.toString(BEARchar)) ){
			return 'y';
			
		}else if(bulgeBranchClass.contains(Character.toString(BEARchar)) ){
			return 'g';
			
		}else if(branchClass.contains(Character.toString(BEARchar)) ){
			return 't';
		}else{
			return '-';
		}
	}
}
