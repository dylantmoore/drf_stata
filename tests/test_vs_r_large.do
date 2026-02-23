*! test_vs_r_large.do
*! Large-scale Stata vs R comparison: n=10K-20K, p=15-30

clear all
set more off
adopath ++ "."

local n_pass = 0
local n_fail = 0
local n_tests = 0

/* ============================================================
 * Scenario A: n=10000, p=20, sparse linear
 * ============================================================ */
display as text _n "=== A: n=10000, p=20, sparse linear ==="
local n_tests = `n_tests' + 1

import delimited using "tests/ref/large_sparse_linear.csv", clear
describe, short

timer on 1
drf y x1-x20, gen(stata_pred) ntrees(300) seed(100)
timer off 1
timer list 1

correlate stata_pred r_pred
local r_vs_r = r(rho)
correlate stata_pred y
local r_vs_y = r(rho)
quietly correlate r_pred y
local r_r_vs_y = r(rho)

display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'
display as text "  Stata vs Y:  r = " as result %7.4f `r_vs_y'
display as text "  R vs Y:      r = " as result %7.4f `r_r_vs_y'

if `r_vs_r' > 0.90 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: Stata-R r = `r_vs_r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Scenario B: n=10000, p=30, nonlinear + interactions
 * ============================================================ */
display as text _n "=== B: n=10000, p=30, nonlinear + interactions ==="
local n_tests = `n_tests' + 1

import delimited using "tests/ref/large_nonlinear_interactions.csv", clear
describe, short

timer on 2
drf y x1-x30, gen(stata_pred) ntrees(300) seed(200)
timer off 2
timer list 2

correlate stata_pred r_pred
local r_vs_r = r(rho)
correlate stata_pred y
local r_vs_y = r(rho)
quietly correlate r_pred y
local r_r_vs_y = r(rho)

display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'
display as text "  Stata vs Y:  r = " as result %7.4f `r_vs_y'
display as text "  R vs Y:      r = " as result %7.4f `r_r_vs_y'

if `r_vs_r' > 0.90 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: Stata-R r = `r_vs_r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Scenario C: n=20000, p=15, heteroskedastic -- mean
 * ============================================================ */
display as text _n "=== C1: n=20000, p=15, heteroskedastic (mean) ==="
local n_tests = `n_tests' + 1

import delimited using "tests/ref/large_heteroskedastic_mean.csv", clear
describe, short

timer on 3
drf y x1-x15, gen(stata_pred) ntrees(300) seed(300)
timer off 3
timer list 3

correlate stata_pred r_pred
local r_vs_r = r(rho)
correlate stata_pred y
local r_vs_y = r(rho)
quietly correlate r_pred y
local r_r_vs_y = r(rho)

display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'
display as text "  Stata vs Y:  r = " as result %7.4f `r_vs_y'
display as text "  R vs Y:      r = " as result %7.4f `r_r_vs_y'

if `r_vs_r' > 0.90 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: Stata-R r = `r_vs_r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Scenario C2: n=20000, p=15, heteroskedastic -- quantile 0.1
 * ============================================================ */
display as text _n "=== C2: n=20000, p=15, heteroskedastic (q10) ==="
local n_tests = `n_tests' + 1

import delimited using "tests/ref/large_heteroskedastic_q10.csv", clear

timer on 4
drf y x1-x15, gen(stata_pred) ntrees(300) seed(300) ///
    functional(quantile) quantile(0.1)
timer off 4
timer list 4

correlate stata_pred r_pred
local r_vs_r = r(rho)
display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'

if `r_vs_r' > 0.90 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: Stata-R r = `r_vs_r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Scenario C3: n=20000, p=15, heteroskedastic -- quantile 0.9
 * ============================================================ */
display as text _n "=== C3: n=20000, p=15, heteroskedastic (q90) ==="
local n_tests = `n_tests' + 1

import delimited using "tests/ref/large_heteroskedastic_q90.csv", clear

timer on 5
drf y x1-x15, gen(stata_pred) ntrees(300) seed(300) ///
    functional(quantile) quantile(0.9)
timer off 5
timer list 5

correlate stata_pred r_pred
local r_vs_r = r(rho)
display as text "  Stata vs R:  r = " as result %7.4f `r_vs_r'

if `r_vs_r' > 0.90 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: Stata-R r = `r_vs_r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================ */
display as text _n "{hline 60}"
if `n_fail' == 0 {
    display as result "Results: `n_pass' / `n_tests' passed (all PASS)"
}
else {
    display as error "Results: `n_pass' / `n_tests' passed, `n_fail' FAILED"
}
display as text "{hline 60}"
