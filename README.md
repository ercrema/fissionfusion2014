# R code for Group Fission-Fusion Agent-Based Model
This repository contains the R scripts for the Fission-Fusion Agent-Based Model explored in the following two pubblications:

Crema, E.R., 2014. [A simulation model of fission-fusion dynamics and long-term settlement change.](https://doi.org/10.1007/s10816-013-9185-4). Journal of Archaeological Method and Theory 21, 385–404.

Crema, E.R., 2015. [Modelling Settlement Rank-Size Fluctuations](https://doi.org/10.1007/978-3-319-00008-4_8), in: Wurzer, G., Kowarik, K., Reschreiter, H. (Eds.), Agent-Based Modeling and Simulation in Archaeology, Advances in Geographic Information Science. Springer International Publishing, pp. 161–181. 

and originally featured on my doctoral thesis:

Crema, E.R. 2013. [Spatial and Temporal Models of Jomon Settlement](https://discovery.ucl.ac.uk/id/eprint/1382589/), Unpublished PhD Thesis, University College London.

where you can find the original code. 

The R scripts contained in this repository has not been updated as the purpose of the repository is to provide direct access to the original code rather than an attempt to improve the model (both in terms of theorethical justifications and computational efficiency). I have however slightly revised the documentation to make this easier to read.

## Setup

To run the model on your machine download the file `fissionfusion.R` on your working directory and type `source('fissionfusion.R')`. There are three versions of the model: 

* `FF()` executes the basic disturbance-free version of the model
* `FF2()` executes the predator-prey version. 
* `FF3()` executes the exogenic disturbance version. 





