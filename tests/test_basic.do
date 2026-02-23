*! test_basic.do -- Basic DRF functionality test
*! Run from the drf_stata directory

clear all
set more off
set seed 42

adopath ++ "."

/* ============================================================
 * Test 1: Basic conditional mean on sysuse auto
 * ============================================================ */
display as text _n "=== Test 1: Basic conditional mean ==="

sysuse auto, clear

drf price mpg weight length, gen(drf_pred) ntrees(100) seed(42)

summarize drf_pred price
correlate drf_pred price

/* ============================================================
 * Test 2: Known linear signal (y = 3*x1 + 2*x2 + noise)
 * ============================================================ */
display as text _n "=== Test 2: Known linear signal ==="

clear
set obs 500
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 3 * x1 + 2 * x2 + rnormal(0, 0.5)

drf y x1 x2 x3, gen(pred_mean) ntrees(200) seed(42)

correlate pred_mean y
local corr_yt = r(rho)

display as text _n "Correlation (pred vs truth): " as result %6.4f `corr_yt'

if `corr_yt' > 0.7 {
    display as result "PASS: OOB conditional mean has reasonable correlation with Y"
}
else {
    display as error "FAIL: OOB correlation with Y = `corr_yt' (expected > 0.7)"
}

/* ============================================================
 * Test 3: Conditional quantile (median)
 * ============================================================ */
display as text _n "=== Test 3: Conditional quantile (median) ==="

drf y x1 x2 x3, gen(pred_q50) ntrees(200) seed(42) ///
    functional(quantile) quantile(0.5) replace

correlate pred_q50 y
local corr_q50 = r(rho)

display as text "Correlation (q50 vs truth): " as result %6.4f `corr_q50'

if `corr_q50' > 0.7 {
    display as result "PASS: Conditional median has reasonable correlation with Y"
}
else {
    display as error "FAIL: Conditional median correlation = `corr_q50'"
}

/* ============================================================
 * Test 4: if/in conditions
 * ============================================================ */
display as text _n "=== Test 4: if/in conditions ==="

drf y x1 x2 x3 if x1 > 0, gen(pred_if) ntrees(100) seed(42)

quietly count if !missing(pred_if) & x1 > 0
local n_pred = r(N)
quietly count if !missing(pred_if) & x1 <= 0
local n_leak = r(N)

display as text "Predictions for x1>0: " as result `n_pred'
display as text "Predictions for x1<=0: " as result `n_leak'

if `n_leak' == 0 & `n_pred' > 0 {
    display as result "PASS: if condition respected"
}
else {
    display as error "FAIL: if condition not respected"
}

/* ============================================================
 * Test 5: replace option
 * ============================================================ */
display as text _n "=== Test 5: replace option ==="

capture drop pred_replace
drf y x1 x2 x3, gen(pred_replace) ntrees(100) seed(42)
drf y x1 x2 x3, gen(pred_replace) ntrees(100) seed(42) replace

display as result "PASS: replace option works"

/* ============================================================
 * Test 6: Nonlinear signal
 * ============================================================ */
display as text _n "=== Test 6: Nonlinear signal (y = x1^2 + sin(x2)) ==="

clear
set obs 1000
set seed 123

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = x1^2 + sin(x2) + rnormal(0, 0.3)

drf y x1 x2 x3, gen(pred_nl) ntrees(300) seed(123)

correlate pred_nl y
local corr_nl = r(rho)

display as text "Correlation (nonlinear pred vs truth): " as result %6.4f `corr_nl'

if `corr_nl' > 0.5 {
    display as result "PASS: DRF captures nonlinear signal"
}
else {
    display as error "FAIL: Nonlinear correlation = `corr_nl' (expected > 0.5)"
}

/* ============================================================
 * Summary
 * ============================================================ */
display as text _n "{hline 50}"
display as text "All basic tests completed."
display as text "{hline 50}"
