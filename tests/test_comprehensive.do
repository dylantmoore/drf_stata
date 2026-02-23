*! test_comprehensive.do -- Comprehensive DRF feature and edge case tests
*! Run from the drf_stata directory

clear all
set more off
adopath ++ "."

local n_pass = 0
local n_fail = 0
local n_tests = 0

/* ============================================================
 * Test 1: Single predictor
 * ============================================================ */
display as text _n "=== Test 1: Single predictor ==="
local n_tests = `n_tests' + 1

clear
set obs 300
set seed 1
gen x = rnormal()
gen y = 2 * x + rnormal(0, 0.5)

drf y x, gen(pred1) ntrees(100) seed(1)

correlate pred1 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.8 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 2: Many predictors (p=10)
 * ============================================================ */
display as text _n "=== Test 2: Many predictors (p=10) ==="
local n_tests = `n_tests' + 1

clear
set obs 500
set seed 2
forvalues j = 1/10 {
    gen x`j' = rnormal()
}
gen y = 3*x1 + 2*x2 - x3 + rnormal(0, 1)

drf y x1-x10, gen(pred2) ntrees(200) seed(2)

correlate pred2 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.7 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 3: nohonesty option
 * ============================================================ */
display as text _n "=== Test 3: nohonesty option ==="
local n_tests = `n_tests' + 1

clear
set obs 300
set seed 3
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 2*x1 + x2 + rnormal(0, 0.5)

drf y x1 x2, gen(pred3) ntrees(100) seed(3) nohonesty

correlate pred3 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.7 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 4: Quantile 0.1 (lower tail)
 * ============================================================ */
display as text _n "=== Test 4: Quantile 0.1 ==="
local n_tests = `n_tests' + 1

clear
set obs 500
set seed 4
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 3*x1 + rnormal(0, 1)

drf y x1 x2, gen(q10) ntrees(200) seed(4) functional(quantile) quantile(0.1)

* q10 should be systematically below the mean
quietly summarize q10
local mean_q10 = r(mean)
quietly summarize y
local mean_y = r(mean)

display "  Mean q10 = " %6.3f `mean_q10' "   Mean y = " %6.3f `mean_y'
if `mean_q10' < `mean_y' {
    display as result "  PASS: q10 mean below y mean"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: q10 not below y mean"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 5: Quantile 0.9 (upper tail)
 * ============================================================ */
display as text _n "=== Test 5: Quantile 0.9 ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(q90) ntrees(200) seed(4) functional(quantile) quantile(0.9)

quietly summarize q90
local mean_q90 = r(mean)

display "  Mean q90 = " %6.3f `mean_q90' "   Mean y = " %6.3f `mean_y'
if `mean_q90' > `mean_y' {
    display as result "  PASS: q90 mean above y mean"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: q90 not above y mean"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 6: Quantile ordering (q10 < q50 < q90)
 * ============================================================ */
display as text _n "=== Test 6: Quantile ordering (q10 < q50 < q90) ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(q50) ntrees(200) seed(4) functional(quantile) quantile(0.5)

quietly summarize q10
local m10 = r(mean)
quietly summarize q50
local m50 = r(mean)
quietly summarize q90
local m90 = r(mean)

display "  Mean q10=" %6.3f `m10' "  q50=" %6.3f `m50' "  q90=" %6.3f `m90'
if `m10' < `m50' & `m50' < `m90' {
    display as result "  PASS: proper quantile ordering"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: quantile ordering violated"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 7: Custom bandwidth
 * ============================================================ */
display as text _n "=== Test 7: Custom bandwidth ==="
local n_tests = `n_tests' + 1

clear
set obs 300
set seed 7
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 2*x1 + rnormal(0, 0.5)

drf y x1 x2, gen(pred7) ntrees(100) seed(7) bandwidth(1.5)

correlate pred7 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 8: Custom mtry
 * ============================================================ */
display as text _n "=== Test 8: Custom mtry ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred8) ntrees(100) seed(7) mtry(1) replace

correlate pred8 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 9: Custom minnodesize
 * ============================================================ */
display as text _n "=== Test 9: Custom minnodesize ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred9) ntrees(100) seed(7) minnodesize(5) replace

correlate pred9 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 10: Custom samplefrac
 * ============================================================ */
display as text _n "=== Test 10: Custom samplefrac ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred10) ntrees(100) seed(7) samplefrac(0.7) replace

correlate pred10 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 11: Custom alpha
 * ============================================================ */
display as text _n "=== Test 11: Custom alpha ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred11) ntrees(100) seed(7) alpha(0.1) replace

correlate pred11 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 12: Custom numfeatures
 * ============================================================ */
display as text _n "=== Test 12: Custom numfeatures ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred12) ntrees(100) seed(7) numfeatures(20) replace

correlate pred12 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 13: Custom honestyfrac
 * ============================================================ */
display as text _n "=== Test 13: Custom honestyfrac ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(pred13) ntrees(100) seed(7) honestyfrac(0.3) replace

correlate pred13 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 14: Small dataset (n=30)
 * ============================================================ */
display as text _n "=== Test 14: Small dataset (n=30) ==="
local n_tests = `n_tests' + 1

clear
set obs 30
set seed 14
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 3*x1 + rnormal(0, 0.5)

drf y x1 x2, gen(pred14) ntrees(50) seed(14) minnodesize(3)

correlate pred14 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.3 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 15: Larger dataset (n=5000)
 * ============================================================ */
display as text _n "=== Test 15: Larger dataset (n=5000) ==="
local n_tests = `n_tests' + 1

clear
set obs 5000
set seed 15
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 2*x1 - x2 + 0.5*x3 + rnormal(0, 1)

drf y x1 x2 x3, gen(pred15) ntrees(200) seed(15)

correlate pred15 y
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.8 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 16: in condition
 * ============================================================ */
display as text _n "=== Test 16: in condition ==="
local n_tests = `n_tests' + 1

drf y x1 x2 x3 in 1/1000, gen(pred16) ntrees(100) seed(16)

quietly count if !missing(pred16) & _n <= 1000
local n_in = r(N)
quietly count if !missing(pred16) & _n > 1000
local n_out = r(N)

display "  Predictions in 1/1000: `n_in'   outside: `n_out'"
if `n_in' > 0 & `n_out' == 0 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 17: Missing values in Y (should be excluded)
 * ============================================================ */
display as text _n "=== Test 17: Missing values ==="
local n_tests = `n_tests' + 1

clear
set obs 200
set seed 17
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 2*x1 + rnormal(0, 0.5)
replace y = . in 1/20

drf y x1 x2, gen(pred17) ntrees(100) seed(17)

quietly count if !missing(pred17)
local n_pred = r(N)
display "  Predictions: `n_pred' out of 200 obs (20 missing y)"
if `n_pred' == 180 {
    display as result "  PASS: missing obs excluded"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: expected 180, got `n_pred'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 18: Missing values in X (should be excluded)
 * ============================================================ */
display as text _n "=== Test 18: Missing X values ==="
local n_tests = `n_tests' + 1

clear
set obs 200
set seed 18
gen x1 = rnormal()
gen x2 = rnormal()
gen y = 2*x1 + rnormal(0, 0.5)
replace x1 = . in 1/10

drf y x1 x2, gen(pred18) ntrees(100) seed(18)

quietly count if !missing(pred18)
local n_pred = r(N)
display "  Predictions: `n_pred' out of 200 obs (10 missing x1)"
if `n_pred' == 190 {
    display as result "  PASS: missing X obs excluded"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: expected 190, got `n_pred'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 19: Reproducibility (same seed = same result)
 * ============================================================ */
display as text _n "=== Test 19: Reproducibility ==="
local n_tests = `n_tests' + 1

clear
set obs 300
set seed 19
gen x1 = rnormal()
gen x2 = rnormal()
gen y = x1 + x2 + rnormal(0, 0.5)

drf y x1 x2, gen(run1) ntrees(100) seed(42)
drf y x1 x2, gen(run2) ntrees(100) seed(42) replace

correlate run1 run2
local r = r(rho)
display "  Correlation between runs: " %8.6f `r'
if `r' > 0.9999 {
    display as result "  PASS: reproducible"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: not reproducible, r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 20: Different seeds = different results
 * ============================================================ */
display as text _n "=== Test 20: Different seeds ==="
local n_tests = `n_tests' + 1

drf y x1 x2, gen(run3) ntrees(100) seed(99) replace

correlate run1 run3
local r = r(rho)
display "  Correlation between different seeds: " %6.4f `r'
if `r' < 0.9999 {
    display as result "  PASS: different seeds give different results"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: seeds not differentiated"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 21: Error on bad functional
 * ============================================================ */
display as text _n "=== Test 21: Bad functional (should error) ==="
local n_tests = `n_tests' + 1

capture drf y x1 x2, gen(bad) ntrees(50) seed(1) functional(variance)
if _rc == 198 {
    display as result "  PASS: bad functional caught (rc=198)"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: expected rc=198, got rc=" _rc
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Test 22: sysuse auto -- a real dataset
 * ============================================================ */
display as text _n "=== Test 22: sysuse auto ==="
local n_tests = `n_tests' + 1

sysuse auto, clear
drf price mpg weight length turn, gen(pred22) ntrees(200) seed(22)

correlate pred22 price
local r = r(rho)
display "  r = " %6.4f `r'
if `r' > 0.5 {
    display as result "  PASS"
    local n_pass = `n_pass' + 1
}
else {
    display as error "  FAIL: r=`r'"
    local n_fail = `n_fail' + 1
}

/* ============================================================
 * Summary
 * ============================================================ */
display as text _n "{hline 50}"
display as result "Results: `n_pass' / `n_tests' passed, `n_fail' failed"
display as text "{hline 50}"
