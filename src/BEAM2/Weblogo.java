package BEAM2;

import java.io.*;

public class Weblogo {
	public void generateLogo(String baseName, int mask, int bestRunIdx,
			String outDad) throws IOException, InterruptedException {

	    String srcname=outDad +"/motifs/" + baseName + "_m" + mask + "_run" + bestRunIdx + "_wl.fa";
	    String outLogo=outDad+"/motifs/" + baseName + "_m" + mask + "_run" + bestRunIdx + "_wl.eps";
	    File tempScript = createTempScript(srcname, outLogo);

	    try {
	        ProcessBuilder pb = new ProcessBuilder("bash", tempScript.toString());
	        //pb.inheritIO();
	        Process process = pb.start();
	        process.waitFor();
	    } finally {
	        tempScript.delete();
	    }
	}

	public File createTempScript(String srcName, String outLogo) throws IOException {
	    File tempScript = File.createTempFile("script", null);
	    String command="./weblogo -a 'ZAQXSWCDEVFRBGTNHY' -f ../../" + srcName + " -o ../../" + outLogo + " -C red ZAQ 'Stem' -C blue XSW 'Loop' -C forestgreen CDE 'InternalLoop' -C orange VFR 'StemBranch' -C DarkOrange B 'Bulge' -C lime G 'BulgeBranch' -C purple T 'Branching' -C limegreen NHY 'InternalLoopBranch'";
	    Writer streamWriter = new OutputStreamWriter(new FileOutputStream(
	            tempScript));
	    PrintWriter printWriter = new PrintWriter(streamWriter);

	    printWriter.println("#!/bin/bash");
	    printWriter.println("cd src/weblogo/");
	    printWriter.println(command);

	    //TODO scrivere script weblogo con srcName
	    printWriter.close();
	    return tempScript;
	}
}
