#' KaH-III Test for Heteroscedasticity (Exact Null Distribution)
#'
#' @description
#' The three-group variance-ratio test for heteroscedasticity in the linear
#' regression model. The data are sorted by a regressor and split into three
#' equal parts; the regression is fitted by ordinary least squares in each part,
#' and the ratio of the largest to the smallest residual mean square is referred
#' to its exact null distribution.
#'
#' @details
#' With \eqn{m = \lfloor n/3 \rfloor} observations per part and \eqn{p}
#' regressors, each residual mean square is
#' \eqn{s_j^2 = \sigma^2 \chi^2_\nu / \nu} with \eqn{\nu = m - p - 1} under
#' homoscedastic normal errors, and the three are independent because the parts
#' share no observations. The statistic
#' \deqn{T = \max_j s_j^2 / \min_j s_j^2}
#' therefore follows Hartley's maximum F-ratio with \eqn{k = 3} groups and
#' \eqn{\nu} degrees of freedom \emph{exactly}, at every sample size (see
#' \code{\link{pfmax}}). This is the three-group extension of the
#' Goldfeld-Quandt test; large values indicate that the error variance changes
#' along the sorting variable.
#'
#' When \eqn{n} is not a multiple of 3, the one or two observations at the part
#' boundaries are left out, exactly as in the original implementation.
#'
#' The test assumes normal errors. With heavy-tailed errors or outliers it
#' rejects too often; use \code{\link{kah.robust.test}} in that case.
#'
#' @param x A fitted \code{lm} object, or a model formula.
#' @param data A data frame, used when \code{x} is a formula.
#' @param order.by The variable the data are sorted by before splitting: the
#'   name of a regressor or of a column of \code{data}, or a numeric vector of
#'   length \eqn{n}. Default: the first regressor.
#'
#' @return An object of class \code{"htest"} with components \code{statistic}
#'   (the ratio \eqn{T}), \code{parameter} (\code{k} and \code{df}),
#'   \code{p.value}, \code{estimate} (the three residual mean squares),
#'   \code{method}, \code{alternative} and \code{data.name}; plus
#'   \code{critical}, the critical values at the 0.5\%, 1\%, 2.5\% and 5\% levels,
#'   and \code{n.used}, the number of observations in the three parts.
#'
#' @references
#' Hartley HO (1950). "The maximum F-ratio as a short-cut test for heterogeneity
#' of variance." \emph{Biometrika}, \bold{37}(3/4), 308-312.
#' \doi{10.1093/biomet/37.3-4.308}
#'
#' Goldfeld SM, Quandt RE (1965). "Some tests for homoscedasticity."
#' \emph{Journal of the American Statistical Association}, \bold{60}(310),
#' 539-547. \doi{10.1080/01621459.1965.10480811}
#'
#' @seealso \code{\link{kah.robust.test}} for the outlier-resistant version,
#'   \code{\link{run.all.het}} to run it next to other tests.
#'
#' @examples
#' # stopping distance of cars: the spread grows with speed
#' kah3.test(dist ~ speed, data = cars)
#'
#' # the same from a fitted model
#' fit <- lm(dist ~ speed, data = cars)
#' kah3.test(fit)$critical
#'
#' @export
kah3.test <- function(x, data = NULL, order.by = NULL) {
  .kah3_core(.het_prepare(x, data, order.by))
}

.kah3_core <- function(d) {
  o <- order(d$ord)
  y <- d$y[o]; Z <- d$Z[o, , drop = FALSE]
  n <- d$n; p <- d$p; m <- floor(n / 3)
  .check_n(m, p, robust = FALSE)

  ms <- vapply(.kah_parts(n), function(i) {
    f <- stats::lm.fit(cbind(1, Z[i, , drop = FALSE]), y[i])
    sum(f$residuals^2) / (length(i) - p - 1)
  }, numeric(1))
  stat <- max(ms) / min(ms)
  nu <- m - p - 1
  levels <- c(0.005, 0.01, 0.025, 0.05)

  structure(list(
    statistic = c(KaH3 = stat),
    parameter = c(k = 3, df = nu),
    p.value = pfmax(stat, nu, 3, lower.tail = FALSE),
    estimate = stats::setNames(ms, c("MSE part 1", "MSE part 2", "MSE part 3")),
    method = "KaH-III test for heteroscedasticity (exact Hartley Fmax null distribution)",
    alternative = paste("error variance changes along", d$ord.name),
    data.name = d$data.name,
    critical = stats::setNames(qfmax(1 - levels, nu, 3), c("0.5%", "1%", "2.5%", "5%")),
    n.used = 3 * m),
    class = "htest")
}
