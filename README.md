# Beam
RNA secondary structure motif discovery.
The standalone version requires the user to provide the input (and eventually background) files in fastB format.
Use the encoder provided on the webserver (beam.uniroma2.it/download).
##The executable version is in .jar format
(if you want to compile the project by yourself, external libraries (commons.io/math) are not included in the repository. Just add them to the project)

### options (in parentheses the default options, if appliable):

-f input file
  A fasta file with additional lines for dotbracket (optional) and BEAR encoding (needed). Informally called fastB (fastBEAR).
 
    >ID
    primary sequence
    dot-bracket
    BEAR notation (encodable with the bear Encoder @ http://beam.uniroma2.it/download)
    >ID2
    ...
    
-g backGround, to compute p-value
  The background file must have the same format as the input, that is a FASTA file with

    >ID
    primary sequence
    dot-bracket
    BEAR notation (encodable with the bear Encoder @ http://beam.uniroma2.it/download)
    >ID2
    ...

-w Min motif Width (10)

-W Max motif Width (50)

-s minSteps (10000)

-S MaxSteps (15000)

-M masks (1) - number of masks(motifs) to be computed 

-R runs (1) - number of runs to be tried before choosing the mask 

**Advanced Options**:

-T Starting Temperature (100)

-r cooling Rate (0.001)

-C clean mode (1,2,3) : 1. Exclude negative partials 2. Exclude partials under the 50% of the partial mean 3. Exclude partials under 90% of the mean. (3)

-b f (false) - weight less (mbr[:][:]=0 ) the branch-branch alignments

-v t Verbose, otherwise mostly silent 

-k keep all the runs, otherwise show best (for R > 1)

-n model_limit (100) – influences the maximum number of structures that form a motif model. It is advisable not to go over this limit, for computational time reasons. The option is here mostly for developers.

-h print help and exit

**Alternative Modes**

--circ treat RNAs as if they were circularized (both sequence and structure, and subopts if appliable)
  you will see a motif that is on the backspliced junction if start and end of a motif position are "switched" (start>end)

## In case you use RNAfold to predict secondary structures
### fold
awk '/^>/ {print; getline; print; getline; print $1}' \<rnafoldOutput\> \> \<encoderReady\>
### encode
java -jar \<encoder.jar\> \<encoderReady\> \<BEAMready\>
then use the \<BEAMready\> file to run BEAM (-f flag)

