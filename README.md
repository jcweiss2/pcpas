## Research code for PCPAS, in JMLR MLHC Conference, 2017

Research code base to run piecewise-constant parametric approximations in survival analysis. Only log logistic parametric distribution is implemented, however the code is modifiable to add others (you provide the parameters, the piecewise-likelihood formulation).

A link to the paper: [piecewise-constant parametric approximations for survival learning](http://www.andrew.cmu.edu/user/jweiss2/2017mlhc_pcpas.pdf)

Synthetic data simulations in 'simulations' folder. Experiments were run on a subpopulation of MIMIC III from [mimic.physionet.org](https://mimic.physionet.org), whose access requires ethics training and and a data use agreement. Details in the paper.


### Running the code
Depedencies include: dplyr, tidyr, compiler, Rcpp, stringr, readr, ggplot2, FAdist. If these packages are not yet installed, use e.g. install.packages("dplyr")
