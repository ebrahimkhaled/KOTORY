#' Effective Degrees of Freedom of the LTS Scale
#'
#' @description
#' The effective degrees of freedom \eqn{\nu^*} used by
#' \code{\link{kah.robust.test}} with \code{method = "fmax"}: the raw least
#' trimmed squares scale of a part with \eqn{m} observations and \eqn{p}
#' regressors behaves approximately like \eqn{\sigma^2\chi^2_{\nu^*}/\nu^*}
#' under normal errors.
#'
#' @details
#' For a scaled chi-square variable, \eqn{Var\{\log(\chi^2_\nu/\nu)\} =
#' \psi'(\nu/2)}, where \eqn{\psi'} is the trigamma function. For every grid
#' point, the LTS scale of one part was simulated under normal errors and
#' \eqn{\nu^*} was obtained by solving \eqn{\psi'(\nu^*/2) =
#' \widehat{Var}(\log \hat\sigma^2_{LTS})}; the same step applied to least
#' squares recovers \eqn{\nu = m - p - 1}. Values between grid points are
#' interpolated on \eqn{\log m}; beyond the largest part size the stable ratio
#' \eqn{\nu^*/\nu} is used. The table is available for \code{alpha} = 0.5, 0.75
#' and 0.9 and \eqn{p = 1, \dots, 5}; for other settings use
#' \code{kah.robust.test(method = "mc")}.
#'
#' @param m Number of observations in each part, \eqn{\lfloor n/3 \rfloor}.
#' @param p Number of regressors (excluding the intercept).
#' @param alpha LTS coverage: 0.5, 0.75 or 0.9.
#'
#' @return The effective degrees of freedom \eqn{\nu^*} (a single number).
#'
#' @examples
#' kah.nu.star(m = 20, p = 1, alpha = 0.75)
#' # ratio to the least squares degrees of freedom
#' kah.nu.star(100, 2, 0.9) / (100 - 2 - 1)
#'
#' @export
kah.nu.star <- function(m, p, alpha = 0.75) {
  tab <- .nu_table
  a <- sort(unique(tab$alpha))
  if (!any(abs(alpha - a) < 1e-8))
    stop("The effective degrees of freedom are tabulated for alpha = ",
         paste(a, collapse = ", "), "; use kah.robust.test(method = \"mc\") for alpha = ",
         alpha, ".", call. = FALSE)
  if (!p %in% tab$p)
    stop("The effective degrees of freedom are tabulated for p = ",
         paste(sort(unique(tab$p)), collapse = ", "),
         " regressors; use kah.robust.test(method = \"mc\") for p = ", p, ".", call. = FALSE)
  t <- tab[abs(tab$alpha - alpha) < 1e-8 & tab$p == p, ]
  t <- t[order(t$m), ]
  if (m < min(t$m))
    stop(sprintf("Parts of %d observations are too small for alpha = %s and p = %d (minimum %d).",
                 m, alpha, p, min(t$m)), call. = FALSE)
  nu <- m - p - 1
  ratio <- if (m > max(t$m)) {
    mean(t$ratio[t$m >= 50])
  } else {
    stats::approx(log(t$m), t$ratio, log(m))$y
  }
  ratio * nu
}
