*********************************************************************************
* Spatial Crime Location Discrete Choice Model 
* Author: Douglas Newball-Ramírez
* Date: 27/06/2022
*********************************************************************************
/* Do-File 1: This Do-file creates descriptive statistics of both datasets */
clear all
set more off
set varabbrev off
graph set window fontface "LM Roman 10"
cd "."

*================================================================================
* 1. Descriptive statistics of the main dataset
*================================================================================
use "crime_barrio_localidad.dta"

* 1.1. Table of simple descriptive statistics
*--------------------------------------------------------------------------------
mat define stats = J(25,5,.)
global crimes crime_v_a_e crime_p_a_e crime_all_a_e
global police tpat_min_e tpat_min
global X meters ss_density near_comercial near_educacion near_parques 			///
	near_policia near_religion near_salud near_servicios_adicionales 			///
	near_transporte closest_camera industry_commerce housing					///
	income_med income_high hs_m treat_hsp

*) Dependent varaiables
local n = 1
foreach y in $crimes {
    local ++n
    sum `y'
	mat stats[`n',1] = r(N)
	mat stats[`n',2] = r(mean)
	mat stats[`n',3] = r(sd)
	mat stats[`n',4] = r(min)
	mat stats[`n',5] = r(max)
}

*) Police presence
local ++n
foreach p in $police {
    local ++n
	sum `p'
	mat stats[`n',1] = r(N)
	mat stats[`n',2] = r(mean)
	mat stats[`n',3] = r(sd)
	mat stats[`n',4] = r(min)
	mat stats[`n',5] = r(max)
    
}


*) Location characteristics	
local ++n
foreach x in $X {
    local ++n
	sum `x'
	mat stats[`n',1] = r(N)
	mat stats[`n',2] = r(mean)
	mat stats[`n',3] = r(sd)
	mat stats[`n',4] = r(min)
	mat stats[`n',5] = r(max)
    
}

frmttable using "../Results/1_summary_statistics", replace tex 					///
	statmat(stats)  sdec(0,2,2,2,2)												///
	rtitles("A. Reported Crimes" \       										///
			"Violent Crimes" \      											///
			"Property Crimes" \      											///
			"Total Crimes" \      												///
			"B. Police Presence" \      										///
			"Avg. police patrol time (minutes)" \      							///
			"Baseline avg. police patrol time (minutes)" \      				///
			"C. Quadrant characteristics" \      								///
			"Avg. longitude of street segments (mt)" \      					///
			"Avg. built squared meters per street segment meter" \      		///
			"Avg. distance to nearest shopping center" \      					///
			"Avg. distance to nearest educational center" \      				///
			"Avg. distance to nearest park/recreational center" \      			///
			"Avg. distance to nearest police station" \      					///
			"Avg. distance to nearest religious/cultural center" \      		///
			"Avg. distance to nearest health center" \      					///
			"Avg. distance to nearest additional services office" \      		///
			"Avg. distance to nearest transport infrastructure" \      			///
			"Avg. distance to closest surveillance camera" \      				///
			"Prop. of street segments zoned for industry/commerce" \      		///
			"Prop. of street segments zoned for housing" \      				///
			"Prop. of middle income street segments" \     						///
			"Prop. of high income street segments" \ 							///
			"Prop. of hot spot stret segments" \								///
			"Prop. of treated hot spot street segments")						///			
	ctitles("", "N", "Mean", "SD", "Min", "Max")					
	
* 1.2. Scatter plot of crime against time
*--------------------------------------------------------------------------------
global titles ""Violent Crimes" "Property Crimes" "Total Crimes""	
tokenize $crimes

local n = 0
foreach y in $titles {
    local ++n
	
	if `n' == 1 {
	    local pos = 4 
	}
	else {
	    local pos = 4.7
	}
	
	
	cap drop log_``n'' 
	cap drop log_patrol
    gen log_``n'' = log(``n'' + 1)
	gen log_patrol = log(tpat_min_e + 1)

	reg log_``n'' log_patrol //, cluster(localidad)
	local b = round(_b[log_patrol],0.01)
	local se = round(_se[log_patrol],0.01)
	
	tw (scatter log_``n'' log_patrol, mlcolor(gs6) mfcolor(gs6*0.5)) 			///
		(lfit log_``n'' log_patrol, lcolor(black) lpattern(_)),					///
		text(`pos' 5.6 "OLS estimate: `b' (`se')") ylabel(, nogrid)				///
		graphregion(color(white)) plotregion(lcolor(black))						///
		xtitle("log(Patrol time (min) + 1)") ytitle("log(`y' + 1)")				///
		legend(order(1 "Data" 2 "Linear fit") c(1)) name(``n'', replace)		///
		title(`y', color(black)) xscale(titlegap(2))
}

grc1leg $crimes, graphregion(color(white)) ring(0) position(4) xsize(4) ysize(1)
gr_edit legend.yoffset = 19
gr_edit legend.xoffset = -10
graph export "../Results/1_scatter_plot.png", replace as(png) width(3900)	
*================================================================================
* End of the Do-file