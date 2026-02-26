{smcl}
{* *! version 0.1.0  26feb2026}{...}
{viewerjumpto "Syntax" "drf##syntax"}{...}
{viewerjumpto "Description" "drf##description"}{...}
{viewerjumpto "Options" "drf##options"}{...}
{viewerjumpto "Examples" "drf##examples"}{...}
{viewerjumpto "Stored results" "drf##results"}{...}
{viewerjumpto "References" "drf##references"}{...}

{title:Title}

{phang}
{bf:drf} {hline 2} Distributional Random Forests for conditional distribution estimation

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:drf}
{depvar} {indepvars}
{ifin}
{cmd:,} {opt gen:erate(newvar)} [{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(newvar)}}name of new variable for predictions{p_end}

{syntab:Functional}
{synopt:{opt func:tional(string)}}functional to compute; {bf:mean} (default) or {bf:quantile}{p_end}
{synopt:{opt q:uantile(#)}}target quantile when {bf:functional(quantile)}; default is {bf:0.5}{p_end}

{syntab:Forest parameters}
{synopt:{opt nt:rees(#)}}number of trees; default is {bf:500}{p_end}
{synopt:{opt seed(#)}}random seed; default is {bf:42}{p_end}
{synopt:{opt mtry(#)}}number of variables tried at each split; default is {bf:0} (auto){p_end}
{synopt:{opt minn:odesize(#)}}minimum node size; default is {bf:15}{p_end}
{synopt:{opt band:width(#)}}kernel bandwidth for weighting; default is {bf:0} (auto){p_end}
{synopt:{opt samplef:rac(#)}}fraction of observations sampled per tree; default is {bf:0.5}{p_end}

{syntab:Honesty}
{synopt:{opt ho:nesty}}use honest splitting (the default){p_end}
{synopt:{opt noho:nesty}}disable honest splitting{p_end}
{synopt:{opt honestyf:rac(#)}}fraction of data for honest estimation; default is {bf:0.5}{p_end}

{syntab:Feature importance}
{synopt:{opt alpha(#)}}significance level for feature importance; default is {bf:0.05}{p_end}
{synopt:{opt numf:eatures(#)}}number of features for importance; default is {bf:10}{p_end}

{syntab:Other}
{synopt:{opt replace}}overwrite {it:newvar} if it already exists{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:drf} fits a distributional random forest to estimate the full conditional
distribution of {depvar} given {indepvars}, then computes a specified functional
(conditional mean or conditional quantile) for each observation.

{pstd}
Distributional random forests extend standard random forests by learning the
entire conditional distribution rather than just the conditional mean.  This
allows prediction of any distributional functional, including quantiles,
variances, and probability intervals.

{pstd}
The method is based on Cevid, Michel, Meinshausen, and Buhlmann (2022).  This
implementation wraps the C++ backend of the R {bf:drf} package by Loris Michel
and Jeffrey Naef ({browse "https://github.com/lorismichel/drf":lorismichel/drf}).

{pstd}
Predictions are stored in the new variable specified by {opt generate()}.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(newvar)} specifies the name of the new variable to store
predictions.  This option is required.

{dlgtab:Functional}

{phang}
{opt functional(string)} specifies the distributional functional to compute.
Options are:

{phang2}{bf:mean} {hline 2} conditional mean of {depvar} given {indepvars}.
This is the default.{p_end}

{phang2}{bf:quantile} {hline 2} conditional quantile of {depvar} given
{indepvars}.  The target quantile is set by {opt quantile()}.{p_end}

{phang}
{opt quantile(#)} specifies the target quantile when {bf:functional(quantile)}
is selected.  The default is {bf:0.5} (median).  Must be between 0 and 1.

{dlgtab:Forest parameters}

{phang}
{opt ntrees(#)} sets the number of trees in the forest.  Default is {bf:500}.
More trees generally improve stability at the cost of computation time.

{phang}
{opt seed(#)} sets the random number seed for reproducibility.
Default is {bf:42}.

{phang}
{opt mtry(#)} sets the number of candidate variables considered at each split.
Default is {bf:0}, which uses the automatic rule (square root of the number of
predictors).

{phang}
{opt minnodesize(#)} sets the minimum number of observations in a terminal
node.  Default is {bf:15}.

{phang}
{opt bandwidth(#)} sets the kernel bandwidth used for weighting.  Default is
{bf:0}, which selects the bandwidth automatically.

{phang}
{opt samplefrac(#)} sets the fraction of observations subsampled for each
tree.  Default is {bf:0.5}.

{dlgtab:Honesty}

{phang}
{opt honesty} enables honest splitting, where separate subsamples are used for
splitting and estimation.  This is the default.

{phang}
{opt nohonesty} disables honest splitting.

{phang}
{opt honestyfrac(#)} sets the fraction of subsampled data reserved for honest
estimation.  Default is {bf:0.5}.  Only relevant when honesty is enabled.

{dlgtab:Feature importance}

{phang}
{opt alpha(#)} sets the significance level for feature importance testing.
Default is {bf:0.05}.

{phang}
{opt numfeatures(#)} sets the number of features to consider for importance.
Default is {bf:10}.

{dlgtab:Other}

{phang}
{opt replace} allows {cmd:drf} to overwrite {it:newvar} if it already exists
in the dataset.

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Conditional mean prediction{p_end}
{phang2}{cmd:. drf price mpg weight length, gen(price_hat)}{p_end}

{pstd}Conditional median prediction{p_end}
{phang2}{cmd:. drf price mpg weight length, gen(price_med) functional(quantile) quantile(0.5)}{p_end}

{pstd}90th percentile prediction{p_end}
{phang2}{cmd:. drf price mpg weight, gen(price_p90) functional(quantile) quantile(0.9)}{p_end}

{pstd}With more trees and a fixed seed{p_end}
{phang2}{cmd:. drf price mpg weight, gen(price_hat2) ntrees(1000) seed(123)}{p_end}

{pstd}Disable honest splitting{p_end}
{phang2}{cmd:. drf price mpg weight, gen(price_hat3) nohonesty}{p_end}

{pstd}Replace an existing variable{p_end}
{phang2}{cmd:. drf price mpg weight, gen(price_hat) replace}{p_end}

{pstd}With if/in restrictions{p_end}
{phang2}{cmd:. drf price mpg weight if foreign == 0, gen(price_dom)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:drf} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations used{p_end}
{synopt:{cmd:r(n_trees)}}number of trees{p_end}
{synopt:{cmd:r(seed)}}random seed{p_end}
{synopt:{cmd:r(mtry)}}number of split candidates (0 = auto){p_end}
{synopt:{cmd:r(min_node)}}minimum node size{p_end}
{synopt:{cmd:r(bandwidth)}}kernel bandwidth (0 = auto){p_end}
{synopt:{cmd:r(alpha)}}significance level for feature importance{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(functional)}}functional used ({bf:mean} or {bf:quantile}){p_end}
{synopt:{cmd:r(depvar)}}dependent variable name{p_end}
{synopt:{cmd:r(indepvars)}}predictor variable names{p_end}
{synopt:{cmd:r(generate)}}name of generated prediction variable{p_end}

{marker references}{...}
{title:References}

{phang}
Cevid, D., L. Michel, J. Naef, N. Meinshausen, and P. Buhlmann. 2022.
Distributional random forests: Heterogeneity adjustment and multivariate
distributional regression.
{it:Journal of Machine Learning Research} 23(333): 1-79.

{phang}
Michel, L. and J. Naef. 2022.
drf: Distributional Random Forests.
R package version 1.1.0.
{browse "https://github.com/lorismichel/drf":https://github.com/lorismichel/drf}

{title:Author}

{pstd}
Dylan Moore, University of Hawaii at Manoa{p_end}
{pstd}
{browse "https://github.com/dylantmoore/drf_stata":https://github.com/dylantmoore/drf_stata}{p_end}
{smcl}
