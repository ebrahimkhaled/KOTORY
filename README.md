# KOTORY

<!-- badges: start -->
[![R-CMD-check](https://github.com/ebrahimkhaled/KOTORY/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ebrahimkhaled/KOTORY/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Robust three-group tests for heteroscedasticity in linear regression.**

The tests sort the data by a regressor, split them into three equal parts, fit
the regression in each part and compare the error scale of the parts.

| Function | What it does |
|---|---|
| `kah3.test()` | Least squares in each part. The ratio of the largest to the smallest residual mean square follows **Hartley's maximum F-ratio with three groups exactly**, so the p-value is exact at every sample size. |
| `kah.robust.test()` | Least trimmed squares (LTS) in each part, so outliers neither create nor hide heteroscedasticity as easily. Three references: an effective-degrees-of-freedom F-max approximation (default), a Monte Carlo reference, or a residual bootstrap for non-normal errors. |
| `run.all.het()` | Runs the KaH tests next to Goldfeld-Quandt, Breusch-Pagan (Koenker), White and the robust modified Goldfeld-Quandt test of Rana, Midi and Imon (2008) and, when **skedastic** is installed, three more recent tests. |
| `pfmax()`, `dfmax()`, `qfmax()`, `rfmax()` | Hartley's maximum F-ratio distribution for any number of groups and any positive degrees of freedom. |

## Installation

```r
# development version
# install.packages("remotes")
remotes::install_github("ebrahimkhaled/KOTORY")
```

## Example

```r
library(KOTORY)

# stopping distance of cars: the spread grows with speed
kah3.test(dist ~ speed, data = cars)
kah.robust.test(dist ~ speed, data = cars)

# all tests side by side
run.all.het(dist ~ speed, data = cars)

# an outlier fools least squares but not the robust test
set.seed(2)
x  <- runif(60)
y0 <- 1 + x + rnorm(60)
y  <- y0; y[which.min(x)] <- 12
c(clean = kah3.test(y0 ~ x)$p.value,       outlier = kah3.test(y ~ x)$p.value)
c(clean = kah.robust.test(y0 ~ x)$p.value, outlier = kah.robust.test(y ~ x)$p.value)
```

## Origin

The three-group tests were proposed in the doctoral thesis of Ahmed El-Kotory
(Alexandria University), where their critical values were tabulated by
simulation. This package replaces those tables: the least squares version
follows Hartley's (1950) maximum F-ratio exactly, and the robust version is
referred to an effective-degrees-of-freedom approximation, a Monte Carlo
reference or a bootstrap.

## References

- Hartley HO (1950). The maximum F-ratio as a short-cut test for heterogeneity of variance. *Biometrika* 37(3/4), 308-312. <https://doi.org/10.1093/biomet/37.3-4.308>
- Goldfeld SM, Quandt RE (1965). Some tests for homoscedasticity. *JASA* 60(310), 539-547. <https://doi.org/10.1080/01621459.1965.10480811>
- Rousseeuw PJ (1984). Least median of squares regression. *JASA* 79(388), 871-880. <https://doi.org/10.1080/01621459.1984.10477105>
- Rana MS, Midi H, Imon AHMR (2008). A robust modification of the Goldfeld-Quandt test for the detection of heteroscedasticity in the presence of outliers. *Journal of Mathematics and Statistics* 4(4), 277-283. <https://doi.org/10.3844/jmssp.2008.277.283>

## Authors

Ahmed El-Kotory and Ebrahim Khaled Ebrahim (maintainer), Department of
Statistics, Faculty of Business, Alexandria University.
