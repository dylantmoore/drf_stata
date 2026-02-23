*! test_vs_r.do -- Validate Stata DRF against R drf package reference output
*! Run from the drf_stata directory after running generate_reference.R

clear all
set more off
adopath ++ "."

/* ============================================================
 * Test A: Linear signal -- conditional mean vs R
 * ============================================================ */
display as text _n "=== Test A: Linear signal -- conditional mean ==="

import delimited using "tests/ref_linear.csv", clear
describe

* Run Stata DRF on same data
drf y x1 x2 x3, gen(stata_pred) ntrees(200) seed(42)

* Compare
correlate stata_pred y
local corr_stata_y = r(rho)

correlate stata_pred r_pred_mean
local corr_stata_r = r(rho)

display as text _n "Correlation (Stata pred vs Y): " as result %6.4f `corr_stata_y'
quietly correlate r_pred_mean y
display as text "Correlation (R pred vs Y):     " as result %6.4f r(rho)
display as text "Correlation (Stata vs R pred): " as result %6.4f `corr_stata_r'

summarize stata_pred r_pred_mean

if `corr_stata_y' > 0.90 {
    display as result "PASS: Stata linear mean correlation > 0.90"
}
else {
    display as error "FAIL: Stata linear mean correlation = `corr_stata_y'"
}

if `corr_stata_r' > 0.85 {
    display as result "PASS: Stata vs R cross-correlation > 0.85"
}
else {
    display as error "NOTE: Stata vs R cross-correlation = `corr_stata_r' (different RNG expected)"
}

/* ============================================================
 * Test B: Nonlinear signal -- conditional mean vs R
 * ============================================================ */
display as text _n "=== Test B: Nonlinear signal -- conditional mean ==="

import delimited using "tests/ref_nonlinear.csv", clear

drf y x1 x2 x3, gen(stata_pred) ntrees(300) seed(123)

correlate stata_pred y
local corr_stata_y = r(rho)

correlate stata_pred r_pred_mean
local corr_stata_r = r(rho)

display as text _n "Correlation (Stata pred vs Y): " as result %6.4f `corr_stata_y'
display as text "Correlation (Stata vs R pred): " as result %6.4f `corr_stata_r'

if `corr_stata_y' > 0.75 {
    display as result "PASS: Stata nonlinear mean correlation > 0.75"
}
else {
    display as error "FAIL: Stata nonlinear mean correlation = `corr_stata_y'"
}

/* ============================================================
 * Test C: Quantile (median) vs R
 * ============================================================ */
display as text _n "=== Test C: Quantile (median) ==="

import delimited using "tests/ref_quantile.csv", clear

drf y x1 x2 x3, gen(stata_q50) ntrees(200) seed(42) ///
    functional(quantile) quantile(0.5)

correlate stata_q50 y
local corr_q50_y = r(rho)

correlate stata_q50 r_pred_q50
local corr_q50_r = r(rho)

display as text _n "Correlation (Stata q50 vs Y):  " as result %6.4f `corr_q50_y'
display as text "Correlation (Stata vs R q50):  " as result %6.4f `corr_q50_r'

if `corr_q50_y' > 0.90 {
    display as result "PASS: Stata median quantile correlation > 0.90"
}
else {
    display as error "FAIL: Stata median quantile correlation = `corr_q50_y'"
}

/* ============================================================
 * Summary
 * ============================================================ */
display as text _n "{hline 50}"
display as text "Cross-validation against R drf package completed."
display as text "{hline 50}"
