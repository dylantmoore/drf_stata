*! test_stress.do
*! Stress test: n=100,000, p=200, complex nonlinear DGP

clear all
set more off
adopath ++ "."

display as text _n "{hline 60}"
display as text "STRESS TEST: n=100,000, p=200"
display as text "{hline 60}"

* ---- Generate data ----
timer on 1
set seed 42
quietly set obs 100000

* Generate 200 predictors
forvalues j = 1/200 {
    quietly gen double x`j' = rnormal()
}

* Complex nonlinear DGP:
*   - 10 linear terms
*   - 5 squared terms
*   - 5 interaction terms
*   - 5 sin/abs terms
*   - heteroskedastic noise
quietly gen double y = ///
    3*x1 - 2*x2 + 1.5*x3 - x4 + 0.8*x5 ///
    + 0.5*x6 - 0.3*x7 + 0.7*x8 - 0.4*x9 + 0.2*x10 ///
    + 2*x11^2 - 1.5*x12^2 + x13^2 - 0.5*x14^2 + 0.3*x15^2 ///
    + x1*x2 - 0.5*x3*x4 + 0.8*x5*x6 - x7*x8 + 0.3*x9*x10 ///
    + sin(3*x16) + abs(x17) - cos(2*x18) + sin(x19)*x20 ///
    + (1 + abs(x1)) * rnormal(0, 0.5)

timer off 1
display as text _n "Data generation:"
timer list 1

describe, short
summarize y, detail

* ---- Run DRF (conditional mean) ----
display as text _n "{hline 60}"
display as text "Running DRF: mean, ntrees=100"
display as text "{hline 60}"

timer on 2
drf y x1-x200, gen(pred_mean) ntrees(100) seed(999)
timer off 2
display as text _n "DRF total time:"
timer list 2

correlate pred_mean y
local r_mean = r(rho)
display as text "  r(pred, y) = " as result %7.4f `r_mean'

* ---- Run DRF (quantile 0.1) ----
display as text _n "{hline 60}"
display as text "Running DRF: quantile(0.1), ntrees=100"
display as text "{hline 60}"

timer on 3
drf y x1-x200, gen(pred_q10) ntrees(100) seed(999) ///
    functional(quantile) quantile(0.1)
timer off 3
display as text _n "DRF quantile(0.1) time:"
timer list 3

correlate pred_q10 y
local r_q10 = r(rho)
display as text "  r(q10_pred, y) = " as result %7.4f `r_q10'

* ---- Summary ----
display as text _n "{hline 60}"
display as text "STRESS TEST RESULTS"
display as text "{hline 60}"
display as text "  n = 100,000   p = 200   trees = 100"
display as text "  Mean prediction r(pred,y) = " as result %7.4f `r_mean'
display as text "  Q10 prediction  r(pred,y) = " as result %7.4f `r_q10'

if `r_mean' > 0.5 & `r_q10' > 0.3 {
    display as result _n "  PASS: Predictions are meaningful at scale"
}
else {
    display as error _n "  WARN: Weak predictions — may need more trees"
}
display as text "{hline 60}"
