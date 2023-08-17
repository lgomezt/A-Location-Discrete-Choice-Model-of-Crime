# A Location Discrete Choice Model of Crime: Police Elasticity and Optimal Deployment
### Alvaro J. Riascos Villegas, Andrés Hoyos, Douglas Newball-Ramírez and Mateo Dulce Rubio
This repository contains all the instructions to replicate the paper ["A Location Discrete Choice Model of Crime: Police Elasticity and Optimal Deployment"](https://github.com/lgomezt/A-Location-Discrete-Choice-Model-of-Crime/blob/main/A%20Location%20Discrete%20Choice%20Modelo%20of%20Crime.pdf) including its results, tables, and figures.

This repository has three main folders. The content of each one is outlined below:
- **Code:** This folder contains four do-files:
  - 1_descriptive_statistics.do: Produces descriptives statistics.
  - 2_homogeneous_model.do: Produces tables and figures related to the estimation of the main model.
  - 3_homogeneous_spatio_temporal.do: Produces tables and figures related to the estimation of an augmented model that includes day/night as a choice.
  - 4_boostrap.do: Produces the Standard Errors neede for statistical inference related to counterfactual scenarios.
- **Data:** This folder contains two dta files:
  - crime_barrio_localidad.dta: Main dataset. Contains information of crimes at the neighborhood leve.
  - crime_barrio_dia_localidad.dta: Secondary dataset. Contains information of crimes at the neighborhood-day/night level. 
- **Results:** Contains all the results produced by the do-files of "Code & Data" folder. The name of each result file is related to each Do-file by the number with which it starts. Namely, 1_... files in this folder are produced by "Code/1_descriptive_statistics.do", 2_... by "Code/2_... .do", and so on and so forth.

## Abstract
<p align = "justify">
Despite the common belief that police presence reduces crime, there is mixed evidence of such causal effects in major Latin America cities. In this work we identify the casual relationship between police presence and criminal events by using a large dataset of a randomized controlled police intervention in Bogotá D.C., Colombia. We use an Instrumental Variables approach to identify the causal effect of interest. Then we consistently estimate a Conditional Logit discrete choice model with aggregate data that allow us to identify agents' utilities for crime location using Two Stage Least Squares. The estimated parameters allow us to compute the police own and cross-elasticities of crime for each of the spatial locations and to evaluate different police patrolling strategies. The elasticity of crime to police presence is, on average across spatial locations, −0.26 for violent crime, −0.38 for property crime and −0.38 for total crime, all statistically significant. Estimates of cross-elasticities are close to zero; however, spillover effects are non-negligible. Counterfactual analysis of different police deployment strategies show, for an optimal allocating algorithm, an average reduction in violent crime of 7.09%, a reduction in property crimes of 8.48% and a reduction in total crimes of 5.15% at no additional cost. These results show the potential efficiency gains of using the model to deploy police resources in the city without increasing the total police time required.
