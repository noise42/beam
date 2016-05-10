# beam
RNA secondary structure motif discovery
##the executable version is in .jar format
external libraries (commons.io/math) are not included in the repository. Just add them to the project

options:

-f input file

-g backGround, to compute p-value

-w Min motif Width (10)

-W Max motif Width (50)

-s minSteps (10000)

-S MaxSteps (15000)

-M masks - number of masks(motifs) to be computed (1)

-R runs - number of runs to be tried before choosing the mask


-T Starting Temperature (100)

-r cooling Rate (0.001)

-o output folder

-C clean mode (1,2,3) : 1. Exclude negative partials 2. Exclude partials under the 50% of the partial mean 3. Exclude partials under 90% 
of the mean. (3)

-b F (false) - weight less (mbr[:][:]=0 ) the branch-branch alignments

-v t Verbose, otherwise mostly silent 

-k keep all the runs, otherwise show best (for R > 1)

-n model_limit (100) – influences the maximum number of structures that form a model. It is advisable not to go over this limit, for computational time reasons

-h print help and exit

