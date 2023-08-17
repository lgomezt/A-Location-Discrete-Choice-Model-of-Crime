# A Location Discrete Choice Model of Crime: Police Elasticity and Optimal Deployment
### Alvaro J. Riascos Villegas, Andrés Hoyos, Douglas Newball-Ramírez and Mateo Dulce Rubio
This repository contains all the instructions to replicate the paper ["A Location Discrete Choice Model of Crime: Police Elasticity and Optimal Deployment"](https://github.com/lgomezt/A-Location-Discrete-Choice-Model-of-Crime/blob/main/A%20Location%20Discrete%20Choice%20Modelo%20of%20Crime.pdf) including its results, tables, and figures.

This repository has three main folders. The content of each one is outlined below:
- **Code:** This folder contains four do-files:
  - 1_descriptive_statistics.do: Produces descriptives statistics.
  - 2_homogeneous_model.do: Produces tables and figures related to the estimation of the main model.
  - 3_homogeneous_spatio_temporal.do: Produces tables and figures related to the estimation of an augmented model that includes day/night as a choice.
  - 4_boostrap.do: Produces the Standard Errors neede for statistical inference related to counterfactual scenarios.
- **Data:**
  - crime_barrio_localidad.dta: Main dataset. Contains information of crimes at the neighborhood leve.
  - crime_barrio_dia_localidad.dta: Secondary dataset. Contains information of crimes at the neighborhood-day/night level. 
- **Results:** Contains all the results produced by the do-files of "Code & Data" folder. The name of each result file is related to each Do-file by the number with which it starts. Namely, 1_... files in this folder are produced by "Code/1_descriptive_statistics.do", 2_... by "Code/2_... .do", and so on and so forth.
