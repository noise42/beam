# beam
RNA secondary structure motif discovery
##the executable version is in .jar format
(if you want to compile the project by yourself, external libraries (commons.io/math) are not included in the repository. Just add them to the project)

options (in parentheses the default options, if appliable):

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


-T Starting Temperature (100)

-r cooling Rate (0.001)

-o output folder (same name as input)

-C clean mode (1,2,3) : 1. Exclude negative partials 2. Exclude partials under the 50% of the partial mean 3. Exclude partials under 90% of the mean. (3)

-b f (false) - weight less (mbr[:][:]=0 ) the branch-branch alignments

-v t Verbose, otherwise mostly silent 

-k keep all the runs, otherwise show best (for R > 1)

-n model_limit (100) – influences the maximum number of structures that form a motif model. It is advisable not to go over this limit, for computational time reasons. The option is here mostly for developers.

-h print help and exit

