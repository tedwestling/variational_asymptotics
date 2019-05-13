# variational_asymptotics
Replication files for analyses and experiments in "Beyond prediction: A framework for inference with variational approximations in mixture models".

To download the NLSY97 data, follow these steps:

1. Download "data/data.NLSY97" from this repository

2. Direct your browser to https://www.nlsinfo.org/investigator/pages/search.jsp?s=NLSY97

3. Under "Upload Tagset (from PC):" click "Choose File" and select data.NLSY97. The "Review Selected Variables" tab should show 53 variables.

4. Click the "Save/Download" tab and then "Advanced Download"

5. Unselect everything except "Tagset", "R Source Code", and "Codebook of Selected Variables". Next to "Filename" type "data".

6. Click download. Save the results to the data/ directory. The files of importance will be called "data.R" and data.dat". 

To run this code, the user will need to have the following R packages installed:
MASS, reshape, ggplot2, plyr, numDeriv, fastGHQuad, mvQuad, lme4, parallel, DescTools, batch, xtable.

Additionally, This code was developed using version 1.1-7 of lme4
         Newer versions of lme4 have changed certain package syntax and warnings,
         and running the code with the current version of lme4 may produce warnings or errors. 
         To install version 1.1-7 of lme4, run the following code, then restart R.
library(devtools)
install_version("lme4", version = "1.1-7", repos = "http://cran.us.r-project.org")

The code files run sequentially, based on the numbering preceding the name of the file. That is, 01_compile_data.R should be run first, followed by 02_run_random_intercepts.R, and so on. Ecah of these files should be sourced in a clean R session. The three simulation files, 05_random_intercept_simulations.R, 07_random_quadratic_simulations.R, and 09_timing_simulations.R, were all run in a cluster environment using the scheduling software slurm.
