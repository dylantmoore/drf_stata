*! drf.ado -- Distributional Random Forests for Stata
*! Version 0.1.0
*! Implements Cevid, Michel, Naf, Meinshausen & Buhlmann (2022)

program define drf, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            FUNCtional(string)                 ///
            Quantile(real 0.5)                 ///
            NTrees(integer 500)                ///
            SEED(integer 42)                   ///
            MTRY(integer 0)                    ///
            MINNodesize(integer 15)            ///
            BANDwidth(real 0)                  ///
            SAMPLEfrac(real 0.5)               ///
            ALPha(real 0.05)                   ///
            NUMFeatures(integer 10)            ///
            HONesty                            ///
            NOHONesty                          ///
            HONestyfrac(real 0.5)              ///
            REPlace                            ///
        ]

    /* ---- Parse functional ---- */
    local func_type 0
    if `"`functional'"' == "" | `"`functional'"' == "mean" {
        local func_type 0
    }
    else if `"`functional'"' == "quantile" {
        local func_type 1
    }
    else {
        display as error `"functional(`functional') not supported; use mean or quantile"'
        exit 198
    }

    /* ---- Parse honesty ---- */
    local do_honesty 1
    if "`nohonesty'" != "" {
        local do_honesty 0
    }
    if "`honesty'" != "" {
        local do_honesty 1
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    /* ---- Parse varlist ---- */
    gettoken depvar indepvars : varlist
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Mark sample ---- */
    marksample touse
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Create output variable ---- */
    quietly gen double `generate' = .

    /* ---- Display header ---- */
    display as text ""
    display as text "Distributional Random Forests"
    display as text "{hline 50}"
    display as text "Dependent variable:  " as result "`depvar'"
    display as text "Predictors:          " as result "`indepvars'"
    display as text "Observations:        " as result `n_use'
    display as text "Trees:               " as result `ntrees'
    if `func_type' == 0 {
        display as text "Functional:          " as result "conditional mean"
    }
    else {
        display as text "Functional:          " as result "quantile(`quantile')"
    }
    display as text "{hline 50}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    cap program drop drf_plugin
    program drf_plugin, plugin using("drf_plugin_`c_os_'.plugin")

    /* ---- Call plugin ---- */
    plugin call drf_plugin `depvar' `indepvars' `generate' ///
        if `touse',                                        ///
        "`ntrees'"                                         ///
        "`seed'"                                           ///
        "`mtry'"                                           ///
        "`minnodesize'"                                    ///
        "`bandwidth'"                                      ///
        "`samplefrac'"                                     ///
        "`alpha'"                                          ///
        "`numfeatures'"                                    ///
        "`func_type'"                                      ///
        "`quantile'"                                       ///
        "`do_honesty'"                                     ///
        "`honestyfrac'"

    /* ---- Store results ---- */
    return scalar N          = `n_use'
    return scalar n_trees    = `ntrees'
    return scalar seed       = `seed'
    return scalar mtry       = `mtry'
    return scalar min_node   = `minnodesize'
    return scalar bandwidth  = `bandwidth'
    return scalar alpha      = `alpha'
    return local  functional   "`functional'"
    return local  depvar       "`depvar'"
    return local  indepvars    "`indepvars'"
    return local  generate     "`generate'"

    /* ---- Summary stats ---- */
    quietly summarize `generate' if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "Predictions written to: " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean:         " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    display as text ""
end
