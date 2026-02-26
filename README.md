# drf_stata

Distributional Random Forests for Stata. Wraps the C++ backend of the R [drf](https://github.com/lorismichel/drf) package by Loris Michel and Jeffrey Naef.

Distributional random forests estimate the full conditional distribution of Y given X, enabling prediction of any distributional functional including conditional means, quantiles, variances, and probability intervals.

## Installation

Platform-specific packages are available via `net install`. Choose the one matching your operating system:

```stata
* macOS:
net install drf_stata_mac, from("https://raw.githubusercontent.com/dylantmoore/drf_stata/main") replace

* Linux:
net install drf_stata_linux, from("https://raw.githubusercontent.com/dylantmoore/drf_stata/main") replace

* Windows:
net install drf_stata_win, from("https://raw.githubusercontent.com/dylantmoore/drf_stata/main") replace

* All platforms (includes all binaries):
net install drf_stata, from("https://raw.githubusercontent.com/dylantmoore/drf_stata/main") replace
```

Requires Stata 14.0 or later.

## Quick Start

```stata
sysuse auto, clear

* Conditional mean prediction
drf price mpg weight length, gen(price_hat)

* Conditional median
drf price mpg weight length, gen(price_med) functional(quantile) quantile(0.5)

* 90th percentile
drf price mpg weight, gen(price_p90) functional(quantile) quantile(0.9)
```

## Help

After installation, type `help drf` in Stata for full documentation.

## References

Cevid, D., L. Michel, J. Naef, N. Meinshausen, and P. Buhlmann. 2022. Distributional random forests: Heterogeneity adjustment and multivariate distributional regression. *Journal of Machine Learning Research* 23(333): 1-79.

## License

GPL-3.0
