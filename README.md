# README #

* hlmer R package
* Version 0.1.0

### Author ###
* Jeffrey A. Dahlke

The hlmer package is a companion to the lme4 and lmerTest packages that automates the formulation of mixed-effects linear equations from user-supplied model characteristics. Output of the hlmer() function is meant to approximate the output of the HLM7 program by Raudenbush, Bryk, and Congdon (2014; Scientific Software International, Inc.) by supplementing lmer() output with statistics such as reliability estimates for random effects, intraclass correlation coefficients (ICCs), and chi-square tests for random-effect variance.

### Installation ###
To install hlmer, make sure you have the devtools library (and its dependencies) installed. Run the following to install those libraries:

install.packages("devtools", dependencies = TRUE)

Once devtools is installed, run the following command to install the hlmer library from GitHub:

devtools::install_github("jadahlke/hlmer")

With hlmer installed, you can load the library:

library(hlmer)

Take hlmer for a test drive! Here's the code to replicate the analysis in Raudenbush and Bryk's (2002) Table 4.5:

hlmer(y_lvl1 = "MATHACH", cluster = "ID", x_lvl1 = "SES",
            x_lvl2 = "SECTOR", center_lvl1 = "cluster", y_lvl2 = "all",
            y_lvl1means = "all", model_type = 4, data = hsb)
