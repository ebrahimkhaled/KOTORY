#' Hartley's Maximum F-Ratio Distribution
#'
#' @description
#' Distribution function, density, quantile function and random generation for
#' Hartley's maximum F-ratio: the ratio of the largest to the smallest of
#' \code{k} independent mean squares, each with \code{df} degrees of freedom,
#' under equal population variances and normal data.
#'
#' @details
#' Let \eqn{S_1, \dots, S_k} be independent \eqn{\chi^2_\nu} variables. The
#' statistic \eqn{F_{max} = \max_j S_j / \min_j S_j} has distribution function
#' \deqn{P(F_{max} \le c) = k \int_0^1 \left[ G\{c\, G^{-1}(u)\} - u \right]^{k-1} du,
#'   \qquad c \ge 1,}
#' where \eqn{G} is the \eqn{\chi^2_\nu} distribution function. The integral is
#' taken on the probability scale, where the integrand is bounded for every
#' \eqn{\nu > 0}, so non-integer degrees of freedom are allowed. The density is
#' obtained by differentiating under the integral sign.
#'
#' For \code{k = 2}, \eqn{P(F_{max} \le c) = P(1/c \le F_{\nu,\nu} \le c)},
#' the two-sided variance-ratio test.
#'
#' The quantile function in \pkg{SuppDists} (\code{qmaxFratio}) is inaccurate in
#' the far upper tail; these functions integrate numerically instead.
#'
#' @param q,x Vector of quantiles (values of the ratio, \eqn{\ge 1}).
#' @param p Vector of probabilities.
#' @param n Number of observations to generate.
#' @param df Degrees of freedom of each mean square (positive, may be non-integer).
#' @param k Number of groups (integer \eqn{\ge 2}). Default 3.
#' @param lower.tail Logical; if \code{TRUE} (default) probabilities are
#'   \eqn{P(F_{max} \le q)}, otherwise \eqn{P(F_{max} > q)}.
#'
#' @return \code{pfmax} gives the distribution function, \code{dfmax} the
#'   density, \code{qfmax} the quantile function and \code{rfmax} random
#'   deviates. Arguments are recycled to a common length.
#'
#' @references
#' Hartley HO (1950). "The maximum F-ratio as a short-cut test for heterogeneity
#' of variance." \emph{Biometrika}, \bold{37}(3/4), 308-312.
#' \doi{10.1093/biomet/37.3-4.308}
#'
#' @examples
#' # Hartley's tabulated 5% and 1% points for k = 3 groups, 4 d.f.: 15.5 and 37
#' qfmax(c(0.95, 0.99), df = 4, k = 3)
#'
#' # p-value of an observed ratio of 6.2 with three groups of 10 d.f.
#' pfmax(6.2, df = 10, lower.tail = FALSE)
#'
#' # k = 2 reduces to the two-sided F test
#' pfmax(3, df = 8, k = 2)
#' pf(3, 8, 8) - pf(1 / 3, 8, 8)
#'
#' @name fmax
NULL

.pfmax1 <- function(c, df, k) {
  if (is.na(c) || is.na(df)) return(NA_real_)
  if (c <= 1) return(0)
  if (is.infinite(c)) return(1)
  g <- function(u) k * pmax(stats::pchisq(c * stats::qchisq(u, df), df) - u, 0)^(k - 1)
  min(1, stats::integrate(g, 0, 1, rel.tol = 1e-10, subdivisions = 2000L)$value)
}

.dfmax1 <- function(c, df, k) {
  if (is.na(c) || is.na(df)) return(NA_real_)
  if (c <= 1 || is.infinite(c)) return(0)
  g <- function(u) {
    x <- stats::qchisq(u, df)
    k * (k - 1) * pmax(stats::pchisq(c * x, df) - u, 0)^(k - 2) * stats::dchisq(c * x, df) * x
  }
  stats::integrate(g, 0, 1, rel.tol = 1e-10, subdivisions = 2000L)$value
}

.check_fmax_args <- function(df, k) {
  if (any(df <= 0, na.rm = TRUE)) stop("'df' must be positive.", call. = FALSE)
  if (any(k < 2 | k != round(k), na.rm = TRUE)) stop("'k' must be an integer >= 2.", call. = FALSE)
}

#' @rdname fmax
#' @export
pfmax <- function(q, df, k = 3, lower.tail = TRUE) {
  .check_fmax_args(df, k)
  out <- mapply(.pfmax1, q, df, k, USE.NAMES = FALSE)
  if (lower.tail) out else 1 - out
}

#' @rdname fmax
#' @export
dfmax <- function(x, df, k = 3) {
  .check_fmax_args(df, k)
  mapply(.dfmax1, x, df, k, USE.NAMES = FALSE)
}

#' @rdname fmax
#' @export
qfmax <- function(p, df, k = 3, lower.tail = TRUE) {
  .check_fmax_args(df, k)
  if (!lower.tail) p <- 1 - p
  mapply(function(pr, v, kk) {
    if (is.na(pr) || is.na(v)) return(NA_real_)
    if (pr < 0 || pr > 1) return(NaN)
    if (pr == 0) return(1)
    if (pr == 1) return(Inf)
    f <- function(lc) .pfmax1(exp(lc), v, kk) - pr
    exp(stats::uniroot(f, c(1e-12, log(1e12)), tol = 1e-12)$root)
  }, p, df, k, USE.NAMES = FALSE)
}

#' @rdname fmax
#' @export
rfmax <- function(n, df, k = 3) {
  .check_fmax_args(df, k)
  df <- rep_len(df, n); k <- rep_len(k, n)
  vapply(seq_len(n), function(i) {
    s <- stats::rchisq(k[i], df[i])
    max(s) / min(s)
  }, numeric(1))
}
