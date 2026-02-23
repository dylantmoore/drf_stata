*! test_vs_r_comprehensive.do
*! Run each R reference scenario through Stata DRF and compare predictions
*! Run from the drf_stata directory after running generate_reference_comprehensive.R

clear all
set more off
adopath ++ "."

local n_pass = 0
local n_fail = 0
local n_tests = 0

program define run_compare
    args name ntrees seed functional qprob opts threshold

    display as text _n "=== `name' ==="

    import delimited using "tests/ref/`name'.csv", clear

    * Count predictors (columns starting with x)
    describe, varlist
    local allvars `r(varlist)'
    local xvars ""
    foreach v of local allvars {
        if substr("`v'", 1, 1) == "x" {
            local xvars `xvars' `v'
        }
    }

    * Build command
    if "`functional'" == "mean" {
        drf y `xvars', gen(stata_pred) ntrees(`ntrees') seed(`seed') `opts'
    }
    else if "`functional'" == "quantile" {
        drf y `xvars', gen(stata_pred) ntrees(`ntrees') seed(`seed') ///
            functional(quantile) quantile(`qprob') `opts'
    }

    * Compare with R
    correlate stata_pred r_pred
    local r_vs_r = r(rho)

    correlate stata_pred y
    local r_vs_y = r(rho)

    quietly correlate r_pred y
    local r_r_vs_y = r(rho)

    display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'
    display as text "  Stata vs Y:  r = " as result %7.4f `r_vs_y'
    display as text "  R vs Y:      r = " as result %7.4f `r_r_vs_y'

    if `r_vs_r' > `threshold' {
        display as result "  PASS (Stata-R correlation > `threshold')"
        c_local test_pass 1
    }
    else {
        display as error "  FAIL (Stata-R correlation = `r_vs_r', threshold = `threshold')"
        c_local test_pass 0
    }
end

* ---- 1: Single predictor, mean ----
local n_tests = `n_tests' + 1
run_compare "single_pred_mean" 100 1 "mean" 0.5 "" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 2: Many predictors (p=10), mean ----
local n_tests = `n_tests' + 1
run_compare "many_pred_mean" 200 2 "mean" 0.5 "" 0.90
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 3: nohonesty ----
local n_tests = `n_tests' + 1
run_compare "nohonesty" 100 3 "mean" 0.5 "nohonesty" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 4: Quantile 0.1 ----
local n_tests = `n_tests' + 1
run_compare "quantile_10" 200 4 "quantile" 0.1 "" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 5: Quantile 0.9 ----
local n_tests = `n_tests' + 1
run_compare "quantile_90" 200 4 "quantile" 0.9 "" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 6: Quantile 0.5 ----
local n_tests = `n_tests' + 1
run_compare "quantile_50" 200 4 "quantile" 0.5 "" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 7: Custom bandwidth ----
local n_tests = `n_tests' + 1
run_compare "custom_bw" 100 7 "mean" 0.5 "bandwidth(1.5)" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 8: Custom minnodesize ----
local n_tests = `n_tests' + 1
run_compare "custom_minnodesize" 100 7 "mean" 0.5 "minnodesize(5)" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 9: Custom samplefrac ----
local n_tests = `n_tests' + 1
run_compare "custom_samplefrac" 100 7 "mean" 0.5 "samplefrac(0.7)" 0.90
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 10: Custom numfeatures ----
local n_tests = `n_tests' + 1
run_compare "custom_numfeatures" 100 7 "mean" 0.5 "numfeatures(20)" 0.95
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 11: Custom honestyfrac ----
local n_tests = `n_tests' + 1
run_compare "custom_honestyfrac" 100 7 "mean" 0.5 "honestyfrac(0.3)" 0.90
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 12: Small dataset (n=30) ----
local n_tests = `n_tests' + 1
run_compare "small_n30" 50 14 "mean" 0.5 "minnodesize(3)" 0.80
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 13: Large dataset (n=5000) ----
local n_tests = `n_tests' + 1
run_compare "large_n5000" 200 15 "mean" 0.5 "" 0.90
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 14: Nonlinear signal ----
local n_tests = `n_tests' + 1
run_compare "nonlinear" 300 123 "mean" 0.5 "" 0.90
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

* ---- 15: Auto-like pattern ----
local n_tests = `n_tests' + 1
run_compare "auto_like" 200 22 "mean" 0.5 "" 0.85
local n_pass = `n_pass' + `test_pass'
local n_fail = `n_fail' + (1 - `test_pass')

/* ============================================================ */
display as text _n "{hline 60}"
if `n_fail' == 0 {
    display as result "Results: `n_pass' / `n_tests' passed (all PASS)"
}
else {
    display as error "Results: `n_pass' / `n_tests' passed, `n_fail' FAILED"
}
display as text "{hline 60}"
