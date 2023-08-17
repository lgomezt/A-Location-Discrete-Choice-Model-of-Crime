*********************************************************************************
* Spatial Crime Location Discrete Choice Model 
* Author: Douglas Newball-Ramírez
* Date: 19/06/2022
*********************************************************************************
/* Do-File 3: This Do-file estimates the spatio-temporal choice model with 
   homogeneous preferences */
clear all
set more off
set varabbrev off
graph set window fontface "LM Roman 10"
cd "."
use "crime_barrio_dia_localidad.dta"
cd "../Results"

*================================================================================
**# 1. Creating relevant variables
*================================================================================
* Market shares
global crimes crime_v_e_ crime_p_e_ crime_all_e_ 
egen pob_total = sum(poblacion)	
foreach x in $crimes {
    gen S_`x' = (`x'+1)/791223.3849
	egen S0_`x' = sum(S_`x')
	replace S0_`x' = 1 - S0_`x'
	sum S_`x'
	sum S0_`x'
}

* Delta terms
encode time, gen(dark)
replace dark = 2 - dark
label drop dark
foreach x in $crimes {
    gen delta_`x' = log(S_`x') - log(S0_`x')
}

* P and Z
gen tpat_dark = tpat_min_e_*dark
gen treat_dark = treat_hsp*dark

*================================================================================
**# 2. Estimating the model
*================================================================================
* 2.1. Selecting covariates: Double Selection
*--------------------------------------------------------------------------------
encode localidad, gen(id_localidad)
global X meters ss_density near_comercial near_educacion near_parques 			///
	near_policia near_religion near_salud near_servicios_adicionales 			///
	near_transporte closest_camera industry_commerce housing					///
	income_med income_high

*) Selecting relevant covariates in terms of P and Z

* Light
lasso linear tpat_min_e_ $X if dark == 0, cluster(localidad)
local a = e(allvars_sel)
global X_ds_day `a'

lasso linear treat_hsp $X if dark == 0, cluster(localidad)
local a = e(allvars_sel)
global X_ds_day $X_ds_day `a'

* Dark	
lasso linear tpat_min_e_ $X if dark == 1, cluster(localidad)
local a = e(allvars_sel)
global X_ds_night `a'

lasso linear treat_hsp $X if dark == 1, cluster(localidad)
local a = e(allvars_sel)
global X_ds_night $X_ds_night `a'

*) Selecting relevant covariates in terms of Crime
* Light
global Y delta_crime_v_e_ delta_crime_p_e_ delta_crime_all_e_
foreach y in $Y {
    lasso linear `y' $X if dark == 0, cluster(localidad)
	local a = e(allvars_sel)
	global X_ds $X_ds_day `a'
}

* Dark
foreach y in $Y {
    lasso linear `y' $X if dark == 1, cluster(localidad)
	local a = e(allvars_sel)
	global X_ds $X_ds_night `a'
}

* 2.2. Estimating the model
*--------------------------------------------------------------------------------
*) Traditional Model
egen dark_localidad = group(dark localidad)
label var tpat_min_e_ "$\alpha$"
global titles ""Violent crimes" "Property crimes" "Total crimes""
global Y crime_v_e_ crime_p_e_ crime_all_e_
tokenize $Y
local n = 0
foreach y in $titles {
    local ++n
	if `n' == 1 local rep_ap = "replace"
	else local rep_ap = "append"

	ivreg2 delta_``n'' (tpat_min_e_ tpat_dark = treat_hsp treat_dark) 			///
		hs_m tpat_min $X_ds_day $X_ds_night i.id_localidad 						///
		dark, cluster(dark_localidad)
	outreg2 using "3_alpha_estimates.tex", keep(dark tpat_min_e_ 				///
		tpat_dark) `rep_ap' tex nocons											///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 			///
		Locality  FE, Yes) label dec(3)
	outreg2 using "3_beta_day_estimates.tex", keep($X_ds_day) `rep_ap' tex 		///
		ctitle(`y') label dec(4) addtext(Controls, Yes, Previos patrol time,	///
		Yes, Locality  FE, Yes)	nocons
	outreg2 using "3_beta_night_estimates.tex", keep($X_ds_night) `rep_ap' tex 	///
		ctitle(`y') label dec(4) addtext(Controls, Yes, Previos patrol time,	///
		Yes, Locality  FE, Yes)	nocons	
}

ivreg2 delta_crime_p_e_ (tpat_min_e_ = treat_hsp)


