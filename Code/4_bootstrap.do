*********************************************************************************
* Spatial Crime Location Discrete Choice Model 
* Author: Douglas Newball-Ramírez
* Date: 21/02/2022
*********************************************************************************
/* Do-File #: This Do-file calculates the uncertainty of several counterfactual
   scenarios */
clear all
set more off
set varabbrev off
graph set window fontface "LM Roman 10"
cd "."
use "crime_barrio_localidad.dta"
cd ".."

*================================================================================
**# 1. Creating relevant variables
*================================================================================
* Market shares
global crimes crime_v_a_e crime_p_a_e crime_all_a_e crime_all_99_u
egen pob_total = sum(poblacion)	
foreach x in $crimes {
    gen S_`x' = (`x'+1)/791223.3849
	egen S0_`x' = sum(S_`x')
	replace S0_`x' = 1 - S0_`x'
	sum S_`x'
	sum S0_`x'
}

* Delta terms
foreach x in $crimes {
    gen delta_`x' = log(S_`x') - log(S0_`x')
}

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
lasso linear tpat_min_e $X, cluster(localidad)
local a = e(allvars_sel)
global X_ds `a'

lasso linear treat_hsp $X, cluster(localidad)
local a = e(allvars_sel)
global X_ds $X_ds `a'

*) Selecting relevant covariates in terms of Crime
global Y delta_crime_v_a_e delta_crime_p_a_e delta_crime_all_a_e
foreach y in $Y {
    lasso linear `y' $X, cluster(localidad)
	local a = e(allvars_sel)
	global X_ds $X_ds `a'
}

* 2.2. Predicted outcomes
*--------------------------------------------------------------------------------
global Y crime_v_a_e crime_p_a_e crime_all_a_e
local n = 0
foreach y in $Y {
	local ++n
	mat `y'_mean = J(6,6,.)
	mat `y'_sum = J(6,6,.)
	
	* Estimating the model	
	ivreg2 delta_`y' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds  				///
		i.id_localidad, cluster(localidad)
	cap drop X_`y'
	predict X_`y', xb
	replace X_`y' =  X_`y' - _b[tpat_min_e]*tpat_min_e
		
	* Observed outcomes	
	sum `y'
	mat `y'_mean[1,1] = r(mean)
	mat `y'_mean[1,2] = r(sd)
	mat `y'_sum[1,1] = r(sum)
	

	* Prediction 1: Base scenario
	tempvar y_hat 
	predict `y_hat'
	replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
	sum `y_hat'
	mat `y'_mean[2,1] = r(mean)
	mat `y'_sum[2,1] = r(sum)
	
	* Prediction 2: Optimal scenario
	preserve
	merge 1:1 cuadrante using resultados.dta
	tempvar y_hat
	gen `y_hat' = _b[tpat_min_e]*sol_`n' + X_`y'
	replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
	sum `y_hat'
	mat `y'_mean[3,1] = r(mean)
	mat `y'_mean[4,1] = `y'_mean[3,1] - `y'_mean[2,1]
	mat `y'_sum[3,1] = r(sum)
	mat `y'_sum[4,1] = `y'_sum[3,1] - `y'_sum[2,1]
	restore
	
	
	* Prediction 3: No Blattman scenario
	tempvar y_hat 
	gen `y_hat' = _b[tpat_min_e]*tpat_min_e_0 + X_`y'
	replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
	sum `y_hat'
	mat `y'_mean[5,1] = r(mean)
	mat `y'_mean[6,1] = `y'_mean[5,1] - `y'_mean[2,1]
	mat `y'_sum[5,1] = r(sum)
	mat `y'_sum[6,1] = `y'_sum[5,1] - `y'_sum[2,1]

}

mat outcomes = crime_v_a_e_mean, crime_p_a_e_mean, crime_all_a_e_mean,			///
	crime_v_a_e_sum, crime_p_a_e_sum, crime_all_a_e_sum

* 2.3. Bootstrap Standard Errors
*-------------------------------------------------------------------------------*
set seed 1
local n = 0
foreach y in $Y {
	local ++n
	mat `y'_mean_boots = J(1000,3,.)
	mat `y'_sum_boots = J(1000,3,.)
	
	forvalues i = 1/1000 {
		preserve		
			merge 1:1 cuadrante using resultados.dta
			bsample _N
			
			* Estimating the model	
			ivreg2 delta_`y' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds  		///
				i.id_localidad, cluster(localidad)
			cap drop X_`y'
			predict X_`y', xb
			replace X_`y' =  X_`y' - _b[tpat_min_e]*tpat_min_e
			
			* Prediction 1: Base scenario
			tempvar y_hat 
			predict `y_hat'
			replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
			sum `y_hat'
			mat `y'_mean_boots[`i',1] = r(mean)
			mat `y'_sum_boots[`i',1] = r(sum)
			
			* Prediction 2: Optimal scenario
			tempvar y_hat
			gen `y_hat' = _b[tpat_min_e]*sol_`n' + X_`y'
			replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
			sum `y_hat'
			mat `y'_mean_boots[`i',2]  = r(mean)
			mat `y'_sum_boots[`i',2]  = r(sum)	
			
			* Prediction 3: No Blattman scenario
			tempvar y_hat 
			gen `y_hat' = _b[tpat_min_e]*tpat_min_e_0 + X_`y'
			replace `y_hat' = exp(`y_hat' + log(S0_`y'))*791223.3849
			sum `y_hat'
			mat `y'_mean_boots[`i',3]  = r(mean)
			mat `y'_sum_boots[`i',3]  = r(sum)
			
		restore
		dis `i'
	}
	
	foreach x in mean sum {
		preserve
			clear
			svmat `y'_`x'_boots
			sum `y'_`x'_boots1
			mat `y'_`x'[2,2] = r(sd)
			sum `y'_`x'_boots2
			mat `y'_`x'[3,2] = r(sd)
			sum `y'_`x'_boots3
			mat `y'_`x'[5,2] = r(sd)
		restore		
	}
}

mat outcomes = crime_v_a_e_mean, crime_p_a_e_mean, crime_all_a_e_mean,			///
	crime_v_a_e_sum, crime_p_a_e_sum, crime_all_a_e_sum
	
*) Hypothesis testing
global titles ""Violent crimes" "Property crimes" "Total crimes""
tokenize $Y
local n = 0
foreach y in $titles {
	local ++n
	foreach x in mean sum {
		preserve
			clear
			svmat ``n''_`x'_boots
			
			*) Hyp. testing 1: Base scenario distribution vs optimal predicted
			local a = ``n''_`x'[3,1]
			local b = ``n''_`x'[5,1]
			
			* Histogram			
			local title = strproper("`x'")
			qui sum ``n''_`x'_boots1
			if "`x'" ==  "mean" {
				local max = round(r(max),0.01)
				local min = round(`a',0.01)
				local range = round((`max' - `min')/10,0.01)
			}
			else {
				local max = round(r(max))
				local min = round(`a')
				local range = round((`max' - `min')/5)
				local format = ", format(%6.0fc)"
			}
			

			hist ``n''_`x'_boots1, xline(`a', lcolor(maroon) lpattern(-))		///
				xline(`b', lcolor(dknavy) lpattern(_))							///
				name(``n''_`x'_1, replace)										///
				graphregion(color(white)) plotregion(lcolor(black))				///
				blcolor(gs8%70) bfcolor(gs11%70) xtitle("`y' (`x')")			///
				xlabel(`min'(`range')`max' `format') xscale(titlegap(3))
			graph export "``n''_`x'_boots1.png", as(png) width(3900) replace
			
			* P-value
			count if ``n''_`x'_boots1 <= `a'	
			mat ``n''_`x'[4,2] = r(N)/1000
			
			count if ``n''_`x'_boots3 >= `b'	
			mat ``n''_`x'[6,2] = r(N)/1000
			
			*) Hyp. testing 2: Distribution comparison
			
			* Histogram
			tw (hist ``n''_`x'_boots2, blcolor(maroon%70) 						///
				bfcolor(maroon*0.7%70))											///
				(hist ``n''_`x'_boots3, blcolor(navy%70) 						///
				bfcolor(navy*0.7%70))											///
				(hist ``n''_`x'_boots1, blcolor(gs8%70) bfcolor(gs11%70)),		///
				name(``n''_`x'_2, replace)										///
				graphregion(color(white)) plotregion(lcolor(black))				///
				xtitle("`y' (`x')")												///
				xlabel(`format') xscale(titlegap(3))							///
				legend(order(3 "Blattman scenario" 1 "Optimal scenario"			///
					2 "No intervention scenario") col(3) symxsize(3))
			graph export "``n''_`x'_boots2.png", as(png) width(3900) replace	
				
			* P-values
			gen n = _n
			reshape long ``n''_`x'_boots, i(n) j(j) 
			
			** Kolmogorov-Smirnov Test
			ksmirnov ``n''_`x'_boots if j != 3, by(j)
			mat ``n''_`x'[4,4] = r(p_2)	
			
			ksmirnov ``n''_`x'_boots if j != 2, by(j)
			mat ``n''_`x'[6,4] = r(p_1)	
			
			** T-test
			ttest ``n''_`x'_boots if j != 3, by(j)
			mat ``n''_`x'[4,5] = r(mu_2) - r(mu_1)
			mat ``n''_`x'[4,6] = r(p_u)
			
			ttest ``n''_`x'_boots if j != 2, by(j)
			mat ``n''_`x'[6,5] = r(mu_2) - r(mu_1)
			mat ``n''_`x'[6,6] = r(p_l)
			
		restore 
	}
}

mat outcomes = crime_v_a_e_mean, crime_v_a_e_sum, crime_p_a_e_mean, 			///
	crime_p_a_e_sum, crime_all_a_e_mean, crime_all_a_e_sum
	
frmttable using "bootstrap.tex", tex statmat(outcomes) substat(1) replace	
*===============================================================================*