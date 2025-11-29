/*******************************************************************************
FULL DO-FILE FOR THESIS EXTENSION:
EM Commodity vs Non-Commodity Extension
Outputs IRFs for: GDP, Exchange Rates, Capital Inflows
*******************************************************************************/

clear all
set more off
set maxvar 8000

* ---------------------------------------------------------
* Load Data
* ---------------------------------------------------------
cd "C:\Users\avanirs\Downloads\495\Data_programs"
use "Data\LP_data.dta", clear

* Drop COVID period
drop if time_q > tq(2019q4)

* Original sample restrictions (same as paper)
drop if IRR_coarse == 1
drop if IRR_coarse == 6

* Fix units for US 12m rate
replace i_treasury_12m_avg_US = i_treasury_12m_avg_US / 100

* LP parameter setup
global lags     = 4
global horizon  = 10
global conf     = 90
global cinorm   = invnormal($conf/100)

* Scale Bauer-Swanson shock by 10 as in paper
replace mps = mps * 10
global instrument "mps"

* ---------------------------------------------------------
* Generate Core Variables
* ---------------------------------------------------------
xtset IFS_code time_q

gen lcpi         = ln(cpi_IFS)
gen inflationIFS = D.lcpi
winsor2 inflationIFS, replace cuts(4 96)

gen ler  = ln(ER_avg_IFS)
gen g_er = D.ler
gen lgdp = ln(rGDPsa_weo)

* Capital flows variable (clean name)
gen cf = all_DB_inflow_AHKS_GDP_2023

* GDP growth & inflation differentials vs US
foreach var of varlist g_gdp_final inflationIFS {
    gen double v1 = `var' if IFS_code == 111
    bys time_q: egen double maxUS = max(v1)
    gen `var'_diff = `var' - maxUS
    drop v1 maxUS
}

xtset IFS_code time_q, quarterly
sort IFS_code time_q

* ---------------------------------------------------------
* Create Lags (RHS)
* ---------------------------------------------------------
local lagvars "lgdp g_er g_gdp_final_diff inflationIFS_diff i_treasury_12m_avg_US mps cf"

foreach var of local lagvars {
    forvalues l = 1/$lags {
        gen `var'_l`l' = L`l'.`var'
    }
}

* ---------------------------------------------------------
* Create Leads for GDP and ER (LHS)
* ---------------------------------------------------------
foreach var of varlist lgdp g_er {
    forvalues h = 0/$horizon {
        gen `var'_`h' = F`h'.`var'
    }
}

* ---------------------------------------------------------
* Create Leads for Capital Flows (LHS)
* ---------------------------------------------------------
forvalues h = 0/$horizon {
    gen cf_`h' = F`h'.cf
}

* ---------------------------------------------------------
* Trim sample to post-1990
* ---------------------------------------------------------
keep if time_q >= tq(1990q1)

* Commodity groups
gen EM_comm    = (EM==1 & CommodityExporters==1)
gen EM_noncomm = (EM==1 & CommodityExporters==0)

* Controls (global + lagged US rates and shocks)
global x "*g_gdp_final_diff_l* *inflationIFS_diff_l* *i_treasury_12m_avg_US_l* *${instrument}_l*"

*******************************************************************************
* LOOP OVER OUTCOMES – CREATE IRF MATRICES
* Outcomes: lgdp, g_er, cf
*******************************************************************************

local outcomes "lgdp g_er cf"
local titles   `"GDP Response" "Exchange Rate Response" "Capital Inflow Response"'

local k = 1
foreach y of local outcomes {

    local this_title : word `k' of `titles'

    if "`y'" == "lgdp" local lhsbase "lgdp"
    if "`y'" == "g_er" local lhsbase "g_er"
    if "`y'" == "cf"   local lhsbase "cf"

    * --------------------------
    * 1. BASE EM (All EMs)
    * --------------------------
    matrix baseEM_`y' = J($horizon+1, 5, .)
    local r = 1

    forvalues h = 0/$horizon {

        if "`y'" == "cf" {
            reghdfe cf_`h' cf_l* $x ${instrument} if EM==1, absorb(IFS_code) vce(cl IFS_code)
        }
        else {
            reghdfe `lhsbase'_`h' `lhsbase'_l* $x ${instrument} if EM==1, absorb(IFS_code) vce(cl IFS_code)
        }

        matrix baseEM_`y'[`r',1] = `h'
        matrix baseEM_`y'[`r',3] = 100*_b[mps]
        matrix baseEM_`y'[`r',4] = 100*(_b[mps] + ${cinorm}*_se[mps])
        matrix baseEM_`y'[`r',5] = 100*(_b[mps] - ${cinorm}*_se[mps])

        local r = `r' + 1
    }

    * --------------------------
    * 2. COMMODITY EMs
    * --------------------------
    matrix EM_comm_mat_`y' = J($horizon+1, 5, .)
    local r = 1

    forvalues h = 0/$horizon {

        if "`y'" == "cf" {
            reghdfe cf_`h' cf_l* $x ${instrument} if EM_comm==1, absorb(IFS_code) vce(cl IFS_code)
        }
        else {
            reghdfe `lhsbase'_`h' `lhsbase'_l* $x ${instrument} if EM_comm==1, absorb(IFS_code) vce(cl IFS_code)
        }

        matrix EM_comm_mat_`y'[`r',1] = `h'
        matrix EM_comm_mat_`y'[`r',3] = 100*_b[mps]
        matrix EM_comm_mat_`y'[`r',4] = 100*(_b[mps] + ${cinorm}*_se[mps])
        matrix EM_comm_mat_`y'[`r',5] = 100*(_b[mps] - ${cinorm}*_se[mps])

        local r = `r' + 1
    }

    * --------------------------
    * 3. NON-COMMODITY EMs
    * --------------------------
    matrix EM_non_mat_`y' = J($horizon+1, 5, .)
    local r = 1

    forvalues h = 0/$horizon {

        if "`y'" == "cf" {
            reghdfe cf_`h' cf_l* $x ${instrument} if EM_noncomm==1, absorb(IFS_code) vce(cl IFS_code)
        }
        else {
            reghdfe `lhsbase'_`h' `lhsbase'_l* $x ${instrument} if EM_noncomm==1, absorb(IFS_code) vce(cl IFS_code)
        }

        matrix EM_non_mat_`y'[`r',1] = `h'
        matrix EM_non_mat_`y'[`r',3] = 100*_b[mps]
        matrix EM_non_mat_`y'[`r',4] = 100*(_b[mps] + ${cinorm}*_se[mps])
        matrix EM_non_mat_`y'[`r',5] = 100*(_b[mps] - ${cinorm}*_se[mps])

        local r = `r' + 1
    }

    local k = `k' + 1
}

*******************************************************************************
* SINGLE-PANEL GRAPHS FOR FINAL PAPER
* (All EM, Commodity EM, Non-Commodity EM)
*******************************************************************************

* -------------------------------
* GDP IRFs
* -------------------------------
clear
svmat baseEM_lgdp, names(b)
rename b1 horizon
rename b3 irf
rename b4 ub
rename b5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(navy) lwidth(medthick) ///
    title("GDP IRF – All Emerging Markets") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in GDP") ///
    yscale(range(-1.5 1)) ///
    ylabel(-1.5(0.25)1)
graph export "GDP_IRF_EM.pdf", replace

clear
svmat EM_comm_mat_lgdp, names(c)
rename c1 horizon
rename c3 irf
rename c4 ub
rename c5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(maroon) lwidth(medthick) ///
    title("GDP IRF – Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in GDP") ///
    yscale(range(-1.5 1)) ///
    ylabel(-1.5(0.25)1)
graph export "GDP_IRF_Commodity.pdf", replace

clear
svmat EM_non_mat_lgdp, names(n)
rename n1 horizon
rename n3 irf
rename n4 ub
rename n5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(blue) lwidth(medthick) ///
    title("GDP IRF – Non-Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in GDP") ///
    yscale(range(-1.5 1)) ///
    ylabel(-1.5(0.25)1)
graph export "GDP_IRF_NonCommodity.pdf", replace

* -------------------------------
* EXCHANGE RATE IRFs
* -------------------------------
clear
svmat baseEM_g_er, names(b)
rename b1 horizon
rename b3 irf
rename b4 ub
rename b5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(navy) lwidth(medthick) ///
    title("Exchange Rate IRF – All Emerging Markets") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in Exchange Rate") ///
    yscale(range(-5 5)) ///
    ylabel(-5(1)5)
graph export "ER_IRF_EM.pdf", replace

clear
svmat EM_comm_mat_g_er, names(c)
rename c1 horizon
rename c3 irf
rename c4 ub
rename c5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(maroon) lwidth(medthick) ///
    title("Exchange Rate IRF – Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in Exchange Rate") ///
    yscale(range(-5 5)) ///
    ylabel(-5(1)5)
graph export "ER_IRF_Commodity.pdf", replace

clear
svmat EM_non_mat_g_er, names(n)
rename n1 horizon
rename n3 irf
rename n4 ub
rename n5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(blue) lwidth(medthick) ///
    title("Exchange Rate IRF – Non-Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent Change in Exchange Rate") ///
    yscale(range(-5 5)) ///
    ylabel(-5(1)5)
graph export "ER_IRF_NonCommodity.pdf", replace

* -------------------------------
* CAPITAL INFLOW IRFs
* -------------------------------
clear
svmat baseEM_cf, names(b)
rename b1 horizon
rename b3 irf
rename b4 ub
rename b5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(navy) lwidth(medthick) ///
    title("Capital Inflow IRF – All Emerging Markets") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent of GDP (Capital Inflows)") ///
    yscale(range(-0.8 0.8)) ///
    ylabel(-0.8(0.2)0.8)
graph export "CF_IRF_EM.pdf", replace

clear
svmat EM_comm_mat_cf, names(c)
rename c1 horizon
rename c3 irf
rename c4 ub
rename c5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(maroon) lwidth(medthick) ///
    title("Capital Inflow IRF – Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent of GDP (Capital Inflows)") ///
    yscale(range(-0.8 0.8)) ///
    ylabel(-0.8(0.2)0.8)
graph export "CF_IRF_Commodity.pdf", replace

clear
svmat EM_non_mat_cf, names(n)
rename n1 horizon
rename n3 irf
rename n4 ub
rename n5 lb
drop if horizon == 0
twoway rarea lb ub horizon, color(gs14) || line irf horizon, lcolor(blue) lwidth(medthick) ///
    title("Capital Inflow IRF – Non-Commodity Exporters") ///
    xtitle("Horizon (Quarters After Shock)") ///
    ytitle("Percent of GDP (Capital Inflows)") ///
    yscale(range(-0.8 0.8)) ///
    ylabel(-0.8(0.2)0.8)
graph export "CF_IRF_NonCommodity.pdf", replace
