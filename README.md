# Railways and Disease Spreading: Applications of Complex Networks
This repository contains the codes developed throughout the final assignment
of the Complex Networks course of the Master in Complex Systems and Biophysics at UB.

## Folders' content

### structural_analysis 

This folder contains the code used to analyse the strcutural properties of a real network [1].
Data can be found in [2], the repository `gb_ptn` must be saved in the same folder as the `structural_analysis` folder. 
In addition, in order to properly draw the map of the Great Britain territory the documents (country_region.dbf, country_region.prj,
country_region.shp, country_region.shx) from [3],[4] must be included in the `structural_analysis` folder. Plots are saved locally in a folder called `images`.

### SIS_model
This folder contain three different codes to study the SIS model and its phase transition. 

#### simulation_new.f90
This code contains the functions to simulate the SIS model. It must be compiled to create the `net2.exe` file.

#### simulation_python.py
This code executes the simulation of the SIS model for the studied parameters (note that not all $\lambda$ values are defined since the
code was modified throughouth the execution of the different simulations.). This code executes the `net2.exe` 
generated with `simulation_new.f90`. Data is saved locally in a folder called `data_results`.

#### plots_and_analysis.py
This code reads the data previously created and performs the analysis of the phase transition. The plots generated are saved
locally in a folder called `plots`.

## Bibliography
[1] R. de Regt, C. von Ferber, Y. Holovatch, and M. Lebovka, *"Public transportation in Great Britain viewed as a complex network."* **Transportmetrica A: Transport Science** **15**, 722–748 (2019).

[2] R. de Regt, *"GB_PTN."* https://bitbucket.org/deregtr/gb_ptn/src/master/ (accessed May. 5, 2026).

[3] Ordnance Survey, *"Boundary-Line."* https://osdatahub.os.uk/data/downloads/open/BoundaryLine (accessed May 8, 2026).

[4] The National Archives, *"Open Government Licence v3.0."* https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/ (accessed Jul. 3, 2026).
