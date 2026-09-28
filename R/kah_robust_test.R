#' Robust KaH Test for Heteroscedasticity (Least Trimmed Squares)
#'
#' @description
#' The outlier-resistant version of \code{\link{kah3.test}}. The data are
#' sorted by a regressor and split into three equal parts; a least trimmed
#' squares (LTS) regression is fitted in each part, and the ratio of the largest
#' to the smallest raw LTS scale is the test statistic. Outliers therefore
#' neither create nor hide heteroscedasticity as easily as with least squares.
#'
#' @details
#' The statistic is \eqn{T = \max_j \hat\sigma_j / \min_j \hat\sigma_j}, a ratio
#' of \emph{scales} (standard deviations), computed with
#' \code{robustbase::ltsReg(alpha = alpha)}. Three references are available:
#' \describe{
#'   \item{\code{"fmax"}}{(default) \eqn{T^2} is referred to Hartley's maximum
#'     F-ratio with three groups and \eqn{\nu^*} effective degrees of freedom.
#'     The LTS scale behaves like a scaled \eqn{\chi^2_{\nu^*}} variable with
#'     \eqn{\nu^* < \nu = m - p - 1}, because trimming loses information;
#'     \eqn{\nu^*} was estimated by simulation for \code{alpha} = 0.5, 0.75 and
#'     0.9, \eqn{p = 1, \dots, 5} and part sizes up to 200 (interpolated between
#'     grid points; the large-\eqn{m} ratio \eqn{\nu^*/\nu} is used beyond). See
#'     \code{\link{kah.nu.star}}.}
#'   \item{\code{"mc"}}{A Monte Carlo reference under homoscedastic normal errors
#'     for the observed regressors: \code{B} data sets are simulated and the
#'     statistic recomputed. Works for any \code{alpha} and any number of
#'     regressors.}
#'   \item{\code{"bootstrap"}}{A residual bootstrap: the standardized LTS
#'     residuals of the three parts are pooled and resampled, which keeps the
#'     error distribution of the data (heavy tails, outliers) under the null
#'     hypothesis. This is the reference to use when the errors may be
#'     non-normal.}
#' }
#' The statistic is invariant to the regression coefficients and to the error
#' scale, so the Monte Carlo and bootstrap samples are generated with zero
#' coefficients.
#'
#' LTS uses random subsampling. For reproducible results the fits are run
#' under \code{seed}; the caller's random-number stream is left unchanged.
#'
#' @inheritParams kah3.test
#' @param alpha LTS coverage: the fraction of each part that is kept (between
#'   0.5 and 1). Lower values resist more outliers but cost power.
#' @param method Reference distribution: \code{"fmax"}, \code{"mc"} or
#'   \code{"bootstrap"}. See Details.
#' @param B Number of Monte Carlo or bootstrap samples.
#' @param seed Integer seed for the LTS subsampling and the resampling, or
#'   \code{NULL} to use the current random-number stream.
#'
#' @return An object of class \code{"htest"} with components \code{statistic},
#'   \code{parameter} (\code{k}, \code{df} and \code{alpha}), \code{p.value},
#'   \code{estimate} (the three LTS scales), \code{method}, \code{alternative},
#'   \code{data.name}, \code{critical} (critical values of \eqn{T} at the 0.5\%,
#'   1\%, 2.5\% and 5\% levels) and \code{n.used}.
#'
#' @references
#' Rousseeuw PJ (1984). "Least median of squares regression." \emph{Journal of
#' the American Statistical Association}, \bold{79}(388), 871-880.
#' \doi{10.1080/01621459.1984.10477105}
#'
#' Rousseeuw PJ, Van Driessen K (2006). "Computing LTS regression for large data
#' sets." \emph{Data Mining and Knowledge Discovery}, \bold{12}(1), 29-45.
#' \doi{10.1007/s10618-005-0024-4}
#'
#' Hartley HO (1950). "The maximum F-ratio as a short-cut test for heterogeneity
#' of variance." \emph{Biometrika}, \bold{37}(3/4), 308-312.
#' \doi{10.1093/biomet/37.3-4.308}
#'
#' @seealso \code{\link{kah3.test}}, \code{\link{run.all.het}},
#'   \code{\link{kah.nu.star}}.
#'
#' @examples
#' kah.robust.test(dist ~ speed, data = cars)
#'
#' # one outlier in a homoscedastic sample
#' set.seed(2)
#' x <- runif(60)
#' y0 <- 1 + x + rnorm(60)
#' y <- y0
#' y[which.min(x)] <- 12
#' c(clean = kah3.test(y0 ~ x)$p.value, outlier = kah3.test(y ~ x)$p.value)
#' # least squares is fooled; the robust statistic does not move
#' c(clean = kah.robust.test(y0 ~ x)$p.value, outlier = kah.robust.test(y ~ x)$p.value)
#'
#' \donttest{
#' # bootstrap reference, safest when the errors may be non-normal
#' kah.robust.test(y ~ x, method = "bootstrap", B = 199)
#' }
#' @export
kah.robust.test <- function(x, data = NULL, order.by = NULL, alpha = 0.75,
                            method = c("fmax", "mc", "bootstrap"), B = 499,
                            seed = 1) {
  method <- match.arg(method)
  .kah_robust_core(.het_prepare(x, data, order.by), alpha, method, B, seed)
}

.kah_robust_core <- function(d, alpha, method, B, seed) {
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0.5 || alpha > 1)
    stop("'alpha' must be a single number between 0.5 and 1.", call. = FALSE)
  o <- order(d$ord)
  y <- d$y[o]; Z <- d$Z[o, , drop = FALSE]
  n <- d$n; p <- d$p; m <- floor(n / 3)
  .check_n(m, p, robust = TRUE)
  levels <- c(0.005, 0.01, 0.025, 0.05)

  .with_seed(seed, {
    L <- .lts3(Z, y, alpha)
    stat <- max(L$scale) / min(L$scale)

    if (method == "fmax") {
      nus <- kah.nu.star(m, p, alpha)
      pval <- pfmax(stat^2, nus, 3, lower.tail = FALSE)
      crit <- sqrt(qfmax(1 - levels, nus, 3))
      ref <- sprintf("Fmax(3, %.2f) reference for the squared ratio", nus)
    } else {
      Zs <- Z
      # a resample whose LTS fit is degenerate in some third (possible with few distinct
      # regressor values and many exactly fitted points) is skipped and counted
      sim <- function(e) tryCatch({ s <- .lts3(Zs, e, alpha)$scale; max(s) / min(s) },
                                  error = function(err) NA_real_)
      Tb <- if (method == "mc") {
        replicate(B, sim(stats::rnorm(n)))
      } else {
        r <- L$resid[is.finite(L$resid)]
        replicate(B, sim(sample(r, n, replace = TRUE)))
      }
      skipped <- sum(!is.finite(Tb))
      Tb <- Tb[is.finite(Tb)]
      if (skipped > 0.1 * B)
        warning(sprintf("%d of %d resamples gave a degenerate LTS fit and were skipped; the reference may be unreliable.",
                        skipped, B), call. = FALSE)
      Bv <- length(Tb)
      pval <- (1 + sum(Tb >= stat)) / (Bv + 1)
      crit <- stats::quantile(Tb, 1 - levels, names = FALSE, type = 8)
      nus <- NA_real_
      ref <- sprintf("%s reference, B = %d%s", if (method == "mc") "Monte Carlo" else "residual bootstrap", Bv,
                     if (skipped) sprintf(" (%d degenerate resamples skipped)", skipped) else "")
    }
  })

  structure(list(
    statistic = c(KaH.robust = stat),
    parameter = if (is.na(nus)) c(k = 3, alpha = alpha) else c(k = 3, df = nus, alpha = alpha),
    p.value = pval,
    estimate = stats::setNames(L$scale, c("LTS scale part 1", "LTS scale part 2", "LTS scale part 3")),
    method = paste0("Robust KaH test for heteroscedasticity (LTS, alpha = ", alpha, "; ", ref, ")"),
    alternative = paste("error variance changes along", d$ord.name),
    data.name = d$data.name,
    critical = stats::setNames(crit, c("0.5%", "1%", "2.5%", "5%")),
    n.used = 3 * m),
    class = "htest")
}
