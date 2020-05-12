# ChukchiPolarBear

This R package houses data, an R script, and a TMB (Template Model Builder) C++ script to recreate analysis of Chukchi Sea aerial survey data.  Once TMB is installed, the analysis can be recreated by installing the package and running the script "./inst/runChukchiPB.R."  Distance sampling models can also be fitted to Russian polar bear detections using the "./inst/DSampling.R" script.  Note that some of the file paths in the scripts may need to be changed to correctly locate the location of the "ChukchiPB_TMB.cpp" file (currently in the inst directory).

Note that the analysis for publication was run using TMB 1.7.15 and R version 3.6.1.

Installation Instructions for TMB
=============
The TMB implementation depends on R >= 3.0.0 and a variety of other tools.

First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Second, please install TMB (Template Model Builder) from https://github.com/kaskr/adcomp

At the moment, TMB can be installed using 

# devtools command to get TMB from GitHub
install_github("kaskr/adcomp/TMB") 
 
It is also available on CRAN.