*********************************************************************************
* Spatial Crime Location Discrete Choice Model 
* Author: Douglas Newball-Ramírez
* Date: 19/06/2022
*********************************************************************************
/* Do-File 2: This Do-file estimates the spatial choice model with homogeneous 
   preferences */
clear all
set more off
set varabbrev off
graph set window fontface "LM Roman 10"
cd "."
use "crime_barrio_localidad.dta"
cd "../Results"
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

* 2.2. Estimating the model
*--------------------------------------------------------------------------------
*) Traditional Model
label var tpat_min_e "$\alpha$"
global titles ""Violent crimes" "Property crimes" "Total crimes""
global Y crime_v_a_e crime_p_a_e crime_all_a_e
tokenize $Y
local n = 0
foreach y in $titles {
    local ++n
	if `n' == 1 local rep_ap = "replace"
	else local rep_ap = "append"
	
	* No controls
	ivreg2 delta_``n'' (tpat_min_e = treat_hsp) hs_m, cluster(localidad)
	outreg2 using "2_alpha_estimates.tex", keep(tpat_min_e) `rep_ap' tex nocons	///
		ctitle(`y') addtext(Controls, No, Previos patrol time, No, 				///
		Locality  FE, No) label dec(3)

	* Double Selected controls
	ivreg2 delta_``n'' (tpat_min_e = treat_hsp) hs_m $X_ds, cluster(localidad)
	outreg2 using "2_alpha_estimates.tex", keep(tpat_min_e) append tex	nocons	///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, No, 			///
		Locality  FE, No) label dec(3)
	
	* Contolling for previous time patrol
	ivreg2 delta_``n'' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds, 			///
		cluster(localidad)
	outreg2 using "2_alpha_estimates.tex", keep(tpat_min_e) append tex nocons	///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 			///
		Locality  FE, No) label dec(3)
	
	* Adding localidad fixed effects
	ivreg2 delta_``n'' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds  			///
		i.id_localidad, cluster(localidad)
	outreg2 using "2_alpha_estimates.tex", keep(tpat_min_e) append tex nocons	///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 			///
		Locality  FE, Yes) label dec(3)
	outreg2 using "2_beta_estimates.tex", keep($X_ds hs_m) `rep_ap' tex nocons	///
		ctitle(`y') label dec(4) addtext(Controls, Yes, Previos patrol time,	///
		Yes, Locality  FE, Yes)	
	gen alpha_``n'' = _b[tpat_min_e]
	predict X_``n'', xb
	replace X_``n'' = X_``n'' -  _b[tpat_min_e]*tpat_min_e
	predict e_``n'', res
}

preserve
keep alpha* X*
export excel using "estimates.xlsx", replace firstrow(var)
restore

* Residual Bootstrap
foreach y in $Y {
	sum e_`y'
	gen c_e_`y' = e_`y'- r(mean)
	mat `y' = J(1050,1,0)
	mat `y'_mean = J(1050,1,0)
	mat `y'_sum = J(1050,1,0)
}


forvalues i = 1/1000 {
	preserve
	bsample 1050
	gen id = _n
	keep id c_e_crime_v_a_e c_e_crime_p_a_e c_e_crime_all_a_e
	rename (c_e_crime_v_a_e c_e_crime_p_a_e c_e_crime_all_a_e)					///
		(c_e_crime_v_a_e_b c_e_crime_p_a_e_b c_e_crime_all_a_e_b)
	tempfile boots
	save `boots', replace
	restore
	
	preserve
	gen id = _n
	merge 1:1 id using `boots'
	foreach y in $Y {
		gen delta_`y'_b = alpha_`y' + X_`y' + c_e_`y'_b
		qui ivreg2 delta_`y'_b (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds  	///
			i.id_localidad
		gen alpha_`y'_`i' = _b[tpat_min_e]
		predict X_`y'_`i', xb	
		replace X_`y'_`i' = X_`y'_`i' - _b[tpat_min_e]*tpat_min_e
		mkmat alpha_`y'_`i' X_`y'_`i'
		mat `y' = `y',alpha_`y'_`i',X_`y'_`i' 
		mat drop alpha_`y'_`i' X_`y'_`i'
		
		predict `y'_hat 
		replace `y'_hat = exp(`y'_hat + log(S0_`y'))*791223.3849
		sum `y'_hat
		mat `y'_mean[`i',1] = r(mean)
		mat `y'_sum[`i',1] = r(sum)
		
	}
	restore
	dis `i'
}

preserve
clear
mat A = crime_p_a_e_mean, crime_v_a_e_mean, crime_all_a_e_mean, crime_p_a_e_sum, ///
	crime_v_a_e_sum, crime_all_a_e_sum
svmat A
drop if _n > 1000
rename (A1 A2 A3 A4 A5 A6) (crime_p_a_e_mean crime_v_a_e_mean					///
	crime_all_a_e_mean crime_p_a_e_sum crime_v_a_e_sum crime_all_a_e_sum)
hist crime_v_a_e_mean, xline(10.345531, lcolor(maroon))
	
	

foreach y in $Y {
	preserve
	clear
	svmat `y'
	drop `y'1
	forvalues i =  1/1000 {
		local j = `i'*2
		rename `y'`j' alpha_`y'_`i'
		local j = `i'*2 + 1
		rename `y'`j' X_`y'_`i'
	}
	tempfile `y'
	gen id = _n
	save ``y'', replace
	restore 
}

preserve
use `crime_v_a_e', clear
merge 1:1 id using `crime_p_a_e'
drop _merge
merge 1:1 id using `crime_all_a_e'
drop _merge
export delimited using "bootstrap.csv", replace
restore

*) Poisson Model
tokenize $Y
local n = 0
foreach y in $titles {
    local ++n
	if `n' == 1 local rep_ap = "replace"
	else local rep_ap = "append"
	
	* GMM
	ivpoisson  gmm ``n'' (tpat_min_e = treat_hsp) hs_m ${X_ds},					///
		vce(cluster id_localidad)
	outreg2 using "2_alpha_estimates_poisson.tex", keep(tpat_min_e) `rep_ap' 	///
		tex nocons ctitle(`y') addtext(Controls, Yes, Previos patrol time, No, 	///
		Locality  FE, No, Method, GMM) label dec(3)
	
	* Control Function
	ivpoisson  cfunc ``n'' (tpat_min_e = treat_hsp) hs_m  ${X_ds},				///
		vce(cluster id_localidad)
	outreg2 using "2_alpha_estimates_poisson.tex", keep(tpat_min_e) append tex	///
		nocons ctitle(`y') addtext(Controls, Yes, Previos patrol time, No, 		///
		Locality  FE, No, Method, CF) label dec(3)
	
	* Control Function + Previous patrol time
	ivpoisson  cfunc ``n'' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds},		///
		vce(cluster id_localidad)
	outreg2 using "2_alpha_estimates_poisson.tex", keep(tpat_min_e) append tex	///
		nocons ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 	///
		Locality  FE, No, Method, CF) label dec(3)
	
	* Control Function + Previous patrol time + Localidad FE
	ivpoisson  cfunc ``n'' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds} 		///
		i.id_localidad, vce(cluster id_localidad)
	outreg2 using "2_alpha_estimates_poisson.tex", keep(tpat_min_e) append tex	///
		nocons ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 	///
		Locality  FE, Yes, Method, CF) label dec(3)
}

*) OLS robustness check
label var tpat_min_e "$\alpha$"
global titles ""Violent crimes" "Property crimes" "Total crimes""
global Y crime_v_a_e crime_p_a_e crime_all_a_e
tokenize $Y
local n = 0
foreach y in $titles {
    local ++n
	if `n' == 1 local rep_ap = "replace"
	else local rep_ap = "append"
	
	* No controls
	reg delta_``n'' tpat_min_e hs_m, cluster(localidad)
	outreg2 using "2_alpha_estimates_ols.tex", keep(tpat_min_e) `rep_ap' tex 	///
		ctitle(`y') addtext(Controls, No, Previos patrol time, No, 				///
		Locality  FE, No) label dec(3) nocons

	* Double Selected controls
	reg delta_``n'' tpat_min_e hs_m $X_ds, cluster(localidad)
	outreg2 using "2_alpha_estimates_ols.tex", keep(tpat_min_e) append tex		///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, No, 			///
		Locality  FE, No) label dec(3) nocons
	
	* Contolling for previous time patrol
	reg delta_``n'' tpat_min_e hs_m tpat_min $X_ds, 							///
		cluster(localidad)
	outreg2 using "2_alpha_estimates_ols.tex", keep(tpat_min_e) append tex 		///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 			///
		Locality  FE, No) label dec(3) nocons
	
	* Adding localidad fixed effects
	reghdfe delta_``n'' tpat_min_e hs_m tpat_min $X_ds, a(id_localidad) 		///
		cluster(localidad)
	outreg2 using "2_alpha_estimates_ols.tex", keep(tpat_min_e) append tex 		///
		ctitle(`y') addtext(Controls, Yes, Previos patrol time, Yes, 			///
		Locality  FE, Yes) label dec(3) nocons

}


* 2.3. Elasticities
*--------------------------------------------------------------------------------
tokenize $Y
local n = 0 
foreach y in $titles {
    local ++n
	if `n' == 1 {
	    local a = "3 -6.5"
		local b = "2.5 -6.5"
		local c = "1.4e5 0.00013"
		local d = "1.2e5 0.00013"
		local e = "0(7.5e4)1.5e5"
	}
	else if `n' == 2 {
	    local a = "2.05 -12.2"
		local b = "1.7 -12.2"
		local c = "7.5e4 0.00021"
		local d = "6.5e4 0.00021"
		local e = "0(4e4)8e4"
	}
	else {
	    local a = "2.3 -12"
		local b = "1.9 -12"
		local c = "7.5e4 0.00021"
		local d = "6.5e4 0.00021"
		local e = "0(4e4)8e4"
	}
	
	*) Own elasticities
    ivreg2 delta_``n'' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds 			///
		i.id_localidad, cluster(localidad)
	cap drop own_e_p own_e_p_se own_e_p_95lb own_e_p_95ub   	
	gen own_e_p = _b[tpat_min_e]*(1-S_crime_p_a_e)*tpat_min_e
	gen own_e_p_se = _se[tpat_min_e]*(1-S_crime_p_a_e)*tpat_min_e
	gen own_e_p_95lb = own_e_p - 1.96*own_e_p_se
	gen own_e_p_95ub = own_e_p + 1.96*own_e_p_se
	mkmat own_e_p own_e_p_95lb own_e_p_95ub own_e_p_se
	mat elasticities = own_e_p\own_e_p
	mat interval = own_e_p_95lb\own_e_p_95ub
	mat se = own_e_p_se\own_e_p_se
	mat e =  elasticities, interval, se 
	preserve
	clear
	svmat e
	sum e1
	local m = r(mean)
	local mr = round(`m',0.001)
	local mr: di %4.3f `mr'
	sum e3
	local s = r(mean)
	local sr = round(`s',0.001)
	local sr: di %4.3f `sr'
	tw (histogram e2, bcolor(gs11)) (histogram e1, bfcolor(gs7) blcolor(gs6)),	///
		graphregion(color(white)) plotregion(lcolor(black)) ylabel(, nogrid)	///
		legend(order(2 "Police elasticity of crime" 1 "95% CI") c(1))	  		///
		xline(`m', lpattern(-) lcolor(black)) text(`a' "Mean = `mr'") 			///
		text(`b' "SE = `sr'") title(`y', color(black)) name(own_e_`n', replace)
	restore
	
	
	*) Cross elasticities
	cap drop cross_e_p cross_e_p_se cross_e_p_95lb cross_e_p_95ub	
	gen cross_e_p = -_b[tpat_min_e]*S_crime_p_a_e*tpat_min_e
	gen cross_e_p_se = _se[tpat_min_e]*S_crime_p_a_e*tpat_min_e
	gen cross_e_p_95lb = cross_e_p - 1.96*cross_e_p_se
	gen cross_e_p_95ub = cross_e_p + 1.96*cross_e_p_se
	mkmat cross_e_p cross_e_p_95lb cross_e_p_95ub cross_e_p_se
	mat elasticities = cross_e_p\cross_e_p
	mat interval = cross_e_p_95lb\cross_e_p_95ub
	mat se = cross_e_p_se\cross_e_p_se
	mat e =  elasticities, interval, se 
	preserve
	clear
	svmat e
	sum e1
	local m = r(mean)
	local mr = round(`m',0.001)
	local mr: di %4.3f `mr'
	sum e3
	local s = r(mean)
	local sr = round(`s',0.001)
	local sr: di %4.3f `sr'
	tw (histogram e2, bcolor(gs11)) (histogram e1, bfcolor(gs7) blcolor(gs6)),	///
		graphregion(color(white)) plotregion(lcolor(black)) ylabel(`e', 		///
		nogrid)	legend(order(2 "Police cross" "elasticity of crime" 1 "95% CI") ///
		c(1)) xline(`m', lpattern(-) lcolor(black)) text(`c' "Mean = `mr'") 	///
		text(`d' "SE = `sr'") title(`y', color(black)) 							///
		name(cross_e_`n', replace) 
	restore
}

grc1leg own_e_1 own_e_2 own_e_3, graphregion(color(white)) ring(0) position(4)
gr_edit legend.yoffset = 17
gr_edit legend.xoffset = -1
graph export "./2_own_elasticities.png", replace as(png) width(3900)

grc1leg cross_e_1 cross_e_2 cross_e_3, graphregion(color(white)) ring(0) 		///
	position(4)
gr_edit legend.yoffset = 17
gr_edit legend.xoffset = -5
graph export "./2_cross_elasticities.png", replace as(png) width(3900)

* 2.4. Policy counterfactual scenarios
*--------------------------------------------------------------------------------
global Y crime_v_a_e crime_p_a_e crime_all_a_e

*) Average crimes
foreach z in $Y {	
	mat `z' = J(15,6,.)
	sum `z'
	mat `z'[1,5] = r(mean)
	mat `z'[1,6] = r(sum)
}	

preserve
clear
import excel "resultados.xlsx", first
gen id_n = _n
save resultados.dta, replace
clear
import excel "resultados_boot.xlsx", first
gen id_n = _n
save resultados_boot.dta, replace
restore

local n = 0
foreach z in $Y {
	local ++n	
	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_`z'))*791223.3849
	sum `y'
	mat `z'[2,5] = r(mean)
	mat `z'[2,6] = r(sum)
	
	preserve
	gen id_n = _n
	merge 1:1 id_n using resultados.dta
	gen y = _b[tpat_min_e]*sol_`n' + X_`z'
	replace y = exp(y + log(S0_`z'))*791223.3849
	sum y
	mat `z'[3,5] = r(mean)	
	mat `z'[3,6] = r(sum)
	restore
	
	preserve
	gen id_n = _n
	merge 1:1 id_n using resultados_boot.dta
	mat `z'_se = J(1000,2,.)
	forvalues i = 1/1000 {
		if "`z'" == "crime_all_a_e" {
			local w = subinstr("`z'","crime_","",.)
		}
		else {
			local w = "`z'"
		}
		
		cap drop y
		gen y = _b[tpat_min_e]*`w'_`i' + X_`z'
		replace y = exp(y + log(S0_`z'))*791223.3849
		sum y
		mat `z'_se[`i',1] = r(mean)
		mat `z'_se[`i',2] = r(sum)		
	}
	restore
	
	preserve 
	clear
	svmat `z'_se
	sum `z'_se1
	mat `z'[4,5] = r(sd)	
	sum `z'_se2
	mat `z'[4,6] = r(sd)	
	restore	
}	

foreach y in $Y {
	preserve
	clear
	svmat `y'_se
	gen id = _n
	tempfile `y'_se
	save ``y'_se.dta', replace
	restore
}

preserve
use `crime_v_a_e_se', replace
merge 1:1 id using `crime_p_a_e_se'
drop _merge
merge 1:1 id using `crime_all_a_e_se'
drop _merge
foreach y in $Y {
	rename (`y'_se1 `y'_se2) (`y'_se_mean `y'_se_sum)
}
save boots_predictions, replace
export excel "boots_predictions.xlsx", first(var)
restore


*) Predictions
foreach z in $Y {	
	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_`z'))*791223.3849
	sum `y'
	mat `z'[2,5] = r(mean)

	mat boots = J(1000,1,.)
	local i = 1
	while `i' <= 1000 {
		dis `i'
		preserve
		bsample 19, cluster(localidad)
		qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}		///
			i.id_localidad, cluster(localidad)
		tempvar y 
		predict `y'
		replace `y' = exp(`y' + log(S0_`z'))*791223.3849
		qui sum `y' if `y' <= 200
		mat boots[`i',1] = r(mean)
		restore
		if boots[`i',1] <= 200 { 
			local ++i // We just count non extreme valued iterations
		}
	}

	preserve
	clear
	svmat boots
	sum boots1
	mat `z'[2,6] = r(sd)
	restore
}

*) Uniform time
foreach z in $Y {
	cap drop total_time
	cap drop total_N
	cap drop time_per_segment
	cap drop tpat_min_unif

	egen total_time = sum(tpat_min_e)
	egen total_N = sum(N)
	gen time_per_segment = total_time/total_N
	gen tpat_min_unif = time_per_segment*N

	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	preserve
	replace tpat_min_e = tpat_min_unif
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_`z'))*791223.3849
	sum `y' if `y' <= 200
	mat `z'[3,5] = r(mean)
	tw kdensity `z' || kdensity `y' if `y' < 200
	restore

	mat boots = J(1000,1,.)
	local i = 1
	while `i' <= 1000 {
		dis `i'
		preserve
		bsample 19, cluster(localidad)
		qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}		///
			i.id_localidad, cluster(localidad)
		replace tpat_min_e = tpat_min_unif	
		tempvar y 
		predict `y'
		replace `y' = exp(`y' + log(S0_`z'))*791223.3849
		qui sum `y' if `y' <= 200
		mat boots[`i',1] = r(mean)
		restore
		if boots[`i',1] <= 200 {
			local ++i // We just count non extreme valued iterations
		}
	}

	preserve
	clear
	svmat boots
	sum boots1
	mat `z'[3,6] = r(sd)
	restore
}

*) Proportional time
foreach z in $Y {
	cap drop n
	cap drop tpat_min_prop

	sort crime_prop_99_u
	gen n = _n
	sort tpat_min_e
	gen tpat_min_prop = tpat_min_e[n]

	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	preserve
	replace tpat_min_e = tpat_min_prop
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_`z'))*791223.3849
	sum `y' if `y' <= 200
	mat `z'[4,5] = r(mean)
	tw kdensity `z' || kdensity `y' if `y' < 200
	restore

	mat boots = J(1000,1,.)
	local i = 1
	while `i' <= 1000 {
		dis `i'
		preserve
		bsample 19, cluster(localidad)
		qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}		///
			i.id_localidad, cluster(localidad)
		replace tpat_min_e = tpat_min_prop
		tempvar y 
		predict `y'
		replace `y' = exp(`y' + log(S0_`z'))*791223.3849
		qui sum `y' if `y' <= 200
		mat boots[`i',1] = r(mean)
		restore
		if boots[`i',1] <= 200 {
			local ++i // We just count non extreme valued iterations
		}
	}

	preserve
	clear
	svmat boots
	sum boots1
	mat `z'[4,6] = r(sd)
	restore
}

*) Reassignment 3
cap drop tmin
cap drop Nmin
cap drop crime_min
cap drop t_diff

egen tmin = min(tpat_min_e)
sort tpat_min_e
gen Nmin = N if _n == 1
ereplace Nmin = mean(Nmin)
gen crime_min = crime_all_99_u if _n == 1
ereplace crime_min = mean(crime_min)

replace tmin = tmin*(N/Nmin)*(crime_all_99_u/crime_min) 	
replace tmin = tpat_min_e if tmin > tpat_min_e					
gen t_diff = tpat_min_e - tmin

foreach z in $Y {
	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	local r = 4
	forvalues j = 1.1(0.1)2.1 {
		local ++r
		sum t_diff
		local diff = r(sum)
		set seed 1
		gsort -crime_all_99_u
		cap drop tpat_min_3
		gen tpat_min_3 = tmin
		local n = 1
		local eval = tpat_min_e[`n']*`j' - tmin[`n']
		while `diff' >= `eval' {
			replace tpat_min_3 = tpat_min_e*`j' if _n == `n'
			local diff = `diff' - `eval' //tpat_min_e[`n']*`j' + tmin[`n']
			local ++n
			dis `n'
			local eval = tpat_min_e[`n']*`j' - tmin[`n']
		}
		local ++n
		replace tpat_min_3 = tpat_min_3 + `diff' if _n == `n' 
		mat `z'[`r',1] = `n'
		mat `z'[`r',3] = (`n'/_N)*100
		 
		preserve
		replace tpat_min_e = tpat_min_3
		tempvar y 
		predict `y'
		replace `y' = exp(`y' + log(S0_`z'))*791223.3849
		sum `y' if `y' <= 200
		mat `z'[`r',5] = r(mean)
		restore
		
		mat boots = J(1000,1,.)
		local i = 1
		while `i' <= 1000 {
			dis `i'
			preserve
			bsample 19, cluster(localidad)
			qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min			///
				${X_ds}	i.id_localidad, cluster(localidad)
			replace tpat_min_e = tpat_min_3
			tempvar y 
			predict `y'
			replace `y' = exp(`y' + log(S0_`z'))*791223.3849
			qui sum `y' if `y' <= 200
			mat boots[`i',1] = r(mean)
			restore
			if boots[`i',1] <= 200 {
				local ++i // We just count non extreme valued iterations
			}
		}

		preserve
		clear
		svmat boots
		sum boots1
		mat `z'[`r',6] = r(sd)
		restore
		
		mat l `z'
	}
}

*) Reassignment 4
foreach z in $Y {
	cap drop tpat_min_4 
	cap drop tpat_min_e_aux 

	gen tpat_min_4 = tmin
	gen tpat_min_e_aux = tpat_min_e	

	qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min $X_ds 			///
		i.id_localidad, cluster(localidad)
	replace tpat_min_e = tpat_min_4

	cap drop c_4 
	cap drop own_e_p_4 
	cap drop rank
	cap drop rank_n

	predict c_4 
	replace c_4 = exp(c_4 + log(S0_`z'))
	gen own_e_p_4 = _b[tpat_min_e]*(1-c_4)*tpat_min_e
	replace c_4 = c_4*791223.3849
	gen rank = c_4 if c_4 <=200	
	
	sum t_diff
	local diff = r(sum)
	
	sum rank
	sort rank
	gen rank_n = _n if rank != .
	sum rank_n
	sum tpat_min_e if rank_n == r(max)
	local a = r(mean)*0.01

	while `diff' >= `a' {	
		sum rank_n
		replace tpat_min_e = tpat_min_e*1.01 if rank_n == r(max)
		local diff = `diff' - `a'		
				
		cap drop c_4 
		cap drop own_e_p_4 
		cap drop rank
		cap drop rank_n
		
		predict c_4 
		replace c_4 = exp(c_4 + log(S0_`z'))
		gen own_e_p_4 = _b[tpat_min_e]*(1-c_4)*tpat_min_e
		replace c_4 = c_4*791223.3849
		
		sum c_4
		gen rank = c_4 if c_4 <= 200
		
		sort rank
		gen rank_n = _n if rank != .
		sum rank_n
		sum tpat_min_e if rank_n == r(max) 
		local a = r(mean)*0.01	
		
		dis `diff'
	}

	cap drop own_e_p_4 
	cap drop rank
	cap drop c_4
	cap drop rank_n
	
	predict c_4 
	replace c_4 = exp(c_4 + log(S0_`z'))
	gen own_e_p_4 = _b[tpat_min_e]*(1-c_4)*tpat_min_e
	replace c_4 = c_4*791223.3849
	gen rank = c_4 if c_4 <= 200	
	
	sort rank
	gen rank_n = _n if rank != .	
	dis `diff'
	sum rank_n
	replace tpat_min_e = tpat_min_e + `diff' if rank_n == r(max) 
	

	replace tpat_min_4 = tpat_min_e
	replace tpat_min_e = tpat_min_e_aux
	cap drop diff
	gen diff = (tpat_min_4 != tmin)
	count if diff == 1
	mat `z'[15,1] == r(N)
	mat `z'[15,3] == (r(N)/_N)*100

	preserve
	replace tpat_min_e = tpat_min_4
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_`z'))*791223.3849
	sum `y' if `y' <= 200
	mat `z'[15,5] = r(mean)
	tw kdensity `z' || kdensity `y' if `y' <= 200
	restore

	mat boots = J(1000,1,.)
	local i = 1
	while `i' <= 1000 {
		dis `i'
		preserve
		bsample 19, cluster(localidad)
		qui ivreg2 delta_`z' (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}		///
			i.id_localidad, cluster(localidad)
		replace tpat_min_e = tpat_min_4
		tempvar y 
		predict `y'
		replace `y' = exp(`y' + log(S0_`z'))*791223.3849
		qui sum `y' if `y' <= 200
		mat boots[`i',1] = r(mean)
		restore
		if boots[`i',1] <= 200 {
			local ++i // We just count non extreme valued iterations
		}
	}

	preserve
	clear
	svmat boots
	sum boots1
	mat `z'[15,6] = r(sd)
	restore
}
/* It might be useful to show the distribution. */

mat C = crime_v_a_e, crime_p_a_e, crime_all_a_e
frmttable using "2_policy_scenario", replace tex 								///
	statmat(C) substat(1) sdec(0,1,2,0,1,2,0,1,2)								///
	rtitles("Observed" \ "" \  													///
			"Predicted" \ "" \ 													///
			"Uniform time" \ "" \ 												///
			"Proportional time" \ "" \ 											///
			"10\% increase" \ "" \ 												///
			"20\% increase" \ "" \ 												///
			"30\% increase" \ "" \ 												///
			"40\% increase" \ "" \ 												///
			"50\% increase" \ "" \ 												///
			"60\% increase" \ "" \ 												///
			"70\% increase" \ "" \ 												///
			"80\% increase" \ "" \ 												///
			"90\% increase" \ "" \ 												///
			"100\% increase" \ "" \ 											///
			"Reassignment 4" \ "" ) 											///
	ctitles("", "Violent Crimes", "", "", 										///
			"Property Crimes", "", "",											///
			"Total Crimes"\														///	
			"","Benefited Q.", "", "Predicted \#",  							///
			"Benefited Q.", "", "Predicted \#",   								///
			"Benefited Q.", "", "Predicted \#" \ 								///
			"", "N", "\%", "Mean (SD) [SE]",  									///
			"N", "\%", "Mean (SD) [SE]", "N", 									///
			"\%", "Mean (SD) [SE]") 											///
	multicol(1,2,3;1,5,3;1,8,3;2,2,2;2,5,2;2,8,2) 												
	
global titles ""Violent Crimes" "Property Crimes" "Total Crimes""	
tokenize $Y
local n = 0
foreach y in $titles {
	local ++n
	if `n' == 1 | `n' == 3 {
		local x = 20
	}
	else {
		local x = 10
	}
		
	preserve
	clear
	svmat ``n''
	drop if _n == 1 | _n == 3 | _n == 4 | _n == 15
	drop ``n''2 ``n''4
	replace ``n''5 = ``n''5
	gen yl = ``n''5 - 1.96*``n''6
	gen yu = ``n''5 + 1.96*``n''6	
	gen N = (_n - 1)*10
	gen V = 10

	tw (rarea yu yl N, color(gs11%70) 											///
		ytitle("Average predicted number of crimes"))							///
		(connected ``n''5 N, msymbol(o) color(black))							///
		(connected ``n''3 N, yaxis(2) msymbol(t) color(maroon) 					///
		ytitle("% of benefited quadrants", axis(2))),							///
		graphregion(color(white)) plotregion(lcolor(black)) ylabel(,nogrid)		///
		xline(`x', lcolor(black) lpattern(-)) 									///
		title(`y', color(black))												///
		xtitle("% increase in police patrol time") 								///
		legend(order(2 "Avg. predicted" "number of crimes" 						///
		1 "95% Confidence" "Interval"											///
		3 "% of benefited" "quadrants") col(1)) 								///
		yscale(titlegap(2) axis(2)) xscale(titlegap(2)) 						///
		name(``n'', replace)
	restore
}	

grc1leg $Y, graphregion(color(white))  ring(0)	position(4) xsize(4) ysize(1)
gr_edit legend.yoffset = 11
gr_edit legend.xoffset = -9
graph export "./2_policy_scenario.png", replace as(png) width(3900)	
		
preserve
clear
svmat crime_v_a_e
svmat crime_p_a_e 
svmat crime_all_a_e
save "./2_matrix_policy.dta", replace
restore	


global dir "G:\Unidades compartidas\Crimen - Elección Discreta\Resultados\2_homogeneous_model"
use "${dir}\2_matrix_policy.dta", clear

global titles ""Violent Crimes" "Property Crimes" "Total Crimes""	
tokenize $Y
local n = 0
foreach y in $titles {
	local ++n
	
	if `n' == 2 {
		local r = "(_n > 5 & _n < 15)"
	}
	else {
		local r = "_n == 5 | (_n > 6 & _n < 15)"
	}
		
	preserve
	keep ``n''*
	drop if _n == 1 | `r'
	sum ``n''5  if _n == 1
	local a = r(mean)
	gen diff = ((``n''5 - `a')/`a')*100
	drop if _n == 1
	gen diff_ub = diff + 1.96*``n''6*100/`a'
	gen diff_lb = diff - 1.96*``n''6*100/`a'
	gen N = _n
	replace N = 5 - N

	tw (rcap diff_ub diff_lb N, color(black) horizontal) 						///
		(scatter N diff, color(black)), 										///
		graphregion(color(white)) plotregion(lcolor(black)) ylabel(,nogrid)		///
		title("`y'", color(black)) xtitle("% change in crimes")					///
		ytitle("") legend(order(2 "% change" 1 "95% CI") col(1))				///
		xline(0, lcolor(maroon) lpattern(-)) ylabel(4 "Uniform" 				///
		3 "Proportional" 2 	"Reassignment 3" 1 "Reasignment 4", angle(0)) 		///
		name(``n'', replace) xscale(titlegap(2))
	restore		
}

grc1leg $Y, graphregion(color(white))  ring(0)	position(4) xsize(4) ysize(1)
gr_edit legend.yoffset = 19
gr_edit legend.xoffset = -8
graph export "${dir}\2_policy_scenario_comparison.png", replace as(png)			///
	width(3900)		
*================================================================================
* End of the Do-file	
	
	
	
/*





	* LATEX
	# delimit ;
	local title1
	"Cuadro 2. Diferencias de medias en el características socioeconómicas";
	local title2 
	"entre niños planeados y no planeados";
	local note1 "***p<0.01, **p<0.05, *p<0.1. Desviación estándar para las"; 
	local note2 "medias y error estándar para la diferencia de medias";
	local note3 "entre paréntesis";
	frmttable using "${outputs}\ed_dif_soci", replace sdec(2) tex        
		statmat(stats) substat(1) annotate(significance) asymbol(*,**,***) 
		title("`title1' `title2'")  
		note("`note1' `note2' `note3'")
		ctitles("","Planeados (N = 475)", "N","No planeados (N = 683)","N",
			"Diferencia (P-NP)")
		rtitles("Sexo masculino" \ "" \
				"Índice de riqueza" \ "" \
				"Personas en el hogar" \ "" \
				"Educación de la madre (años)" \ "" \
				"Edad de la madre (años)" \ "" \
				"Vocabulario materno (TVIP)" \ "" \
				"Depresión CESD10" \ "" \
				"Extraversión (Puntaje Z)" \ "" \
				"Agradabilidad (Puntaje Z)" \ "" \
				"Conciencia (Puntaje Z)" \ "" \
				"Neuroticismo (Puntaje Z)" \ "" \
				"Mente abierta (Puntaje Z)" \ "" \
				"Autoeficacia" \ "" \
				"Apoyo social (Puntaje Z)" \ "" \
				"Asignación al tratamiento" \ "" \
				"Proporción de embarazos pasados perdidos" \ "");
	# delimit cr
	












cap program drop predictse
program define predictse, eclass
	qui ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_crime_p_a_e))*791223.3849
	sum `y'
end
set seed 1234
bootstrap r(mean), reps(1000): predictse
mat crime_p_a_e[2,1] = e(b)[1,1]
mat crime_p_a_e[2,2] = e(se)[1,1]


cap program drop predictse
program define predictse, eclass
	qui ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
	preserve
	replace tpat_min_e = tpat_min_unif
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_crime_p_a_e))*791223.3849
	sum `y'
	restore
end
set seed 1234
bootstrap r(mean), reps(1000): predictse
mat crime_p_a_e[2,1] = e(b)[1,1]
mat crime_p_a_e[2,2] = e(se)[1,1]
drop tpat_min_unif




ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
	 i.id_localidad, cluster(localidad)
preserve


preserve
qui ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
		i.id_localidad, cluster(localidad)
replace tpat_min_e = tpat_min_unif
sum tpat_min_e
tempvar y yse
predict `y'
replace `y' = exp(`y' + log(S0_crime_p_a_e))*791223.3849
sum `y'
mat crime_p_a_e[3,1] = r(mean)
predict `yse', stdp
replace `yse' = `y'*`yse'
sum `yse'
mat crime_p_a_e[3,2] = r(sd)
restore

sum tpat_min_e
tempvar y yse
predict `y'
replace `y' = exp(`y' + log(S0_crime_p_a_e))*791223.3849
sum `y'
mat crime_p_a_e[3,1] = r(mean)
predict `yse', stdp
replace `yse' = `y'*`yse'
sum `yse'
mat crime_p_a_e[3,2] = r(sd)
restore
















cap program drop predictse
program define predictse, eclass
	qui ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min ${X_ds}			///
	 i.id_localidad, cluster(localidad)
	tempvar y 
	predict `y'
	replace `y' = exp(`y' + log(S0_crime_p_a_e))*791223.3849
	sum `y'
end
set seed 1234
bootstrap r(mean), reps(1000): predictse




global Y crime_v_a_e crime_p_a_e crime_all_a_e
local n = 0
foreach y in $Y {
    local ++n
	if `n' == 1 local rep_ap = "replace"
	else local rep_ap = "append"
	
	* Traditional Model
	ivreg2 delta_`y' (tpat_min_e = treat_hsp) hs_m $X_ds, cluster(localidad)
	outreg2 using "2_alpha_estimates.tex", keep(tpat_min_e) `repap'
	
	* Poisson Model
    ivpoisson gmm `y' (tpat_min_e = treat_hsp) hs_m $X_ds,				 		///
		vce(cluster id_localidad)
}




















encode time, gen(dark)
replace dark = 2 - dark
label drop dark

ivreg2 delta_crime_p_a_e (tpat_min_e_ c.tpat_min_e_#dark = treat_hsp c.treat_hsp#dark) hs_m $x1, cluster(cuadrante) first


ivreg2 delta_crime_v_a_e (tpat_min_e = treat_hsp) crime_all_99_u tpat_min, cluster(cuadrante)
ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m $x1, cluster(cuadrante)
ivreg2 delta_crime_all_a_e (tpat_min_e = treat_hsp) hs_m $x1, cluster(cuadrante)








ivreg2 crime_p_a_e (tpat_min_e = treat_hsp) hs_m tpat_min crime_all_99_u
encode barrio, gen(id_barrio)
ivreg2 delta_crime_p_a_e (tpat_min_e = treat_hsp) hs_m 
ivpoisson gmm crime_p_e_ (tpat_min_e_ = treat_hsp) hs_m