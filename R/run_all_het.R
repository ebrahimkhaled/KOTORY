#' Run the KaH Tests Next to Other Heteroscedasticity Tests
#'
#' @description
#' Runs the proposed tests and a set of established heteroscedasticity tests on
#' the same model in one call, and returns one row per test. The established
#' tests are provided for comparison and are attributed to their authors; they
#' are not claimed as original to this package.
#'
#' @details
#' Tests always run:
#' \describe{
#'   \item{KaH-III}{\code{\link{kah3.test}}, exact Hartley reference.}
#'   \item{KaH robust}{\code{\link{kah.robust.test}} with \code{alpha} and the
#'     \code{"fmax"} reference.}
#'   \item{Goldfeld-Quandt}{Goldfeld and Quandt (1965): data sorted by
#'     \code{order.by}, the middle third omitted, ratio of the residual mean
#'     squares of the outer thirds, two-sided F reference.}
#'   \item{Breusch-Pagan (Koenker)}{Koenker's (1981) studentized version of the
#'     Breusch and Pagan (1979) test: \eqn{nR^2} of the squared residuals on the
#'     regressors, \eqn{\chi^2_p} reference.}
#'   \item{White}{White (1980): \eqn{nR^2} of the squared residuals on the
#'     regressors, their squares and cross products.}
#'   \item{MGQ (Rana et al.)}{The robust modified Goldfeld-Quandt test of Rana,
#'     Midi and Imon (2008): median squared deletion residuals of the outer
#'     thirds after an LTS outlier screen, F reference as proposed by the
#'     authors (two-sided).}
#' }
#' Optional tests (\code{include.optional = TRUE}, need \pkg{skedastic}):
#' Wilcox and Keselman (2006), Zhou, Song and Thompson (2015) and Li and Yao
#' (2019). Optional: \code{bootstrap = TRUE} adds the robust KaH test with the
#' residual-bootstrap reference.
#'
#' One failing test never stops the run; its row gets \code{NA} and a note.
#'
#' @inheritParams kah.robust.test
#' @param include.optional Run the tests from \pkg{skedastic} if it is
#'   installed.
#' @param bootstrap Also run the robust KaH test with the bootstrap reference
#'   (slower).
#'
#' @return A data frame of class \code{"het_battery"} with columns
#'   \code{Test}, \code{Group}, \code{Statistic}, \code{df}, \code{p_value} and
#'   \code{Note}. It prints grouped, with the p-values formatted.
#'
#' @references
#' Breusch TS, Pagan AR (1979). "A simple test for heteroscedasticity and
#' random coefficient variation." \emph{Econometrica}, \bold{47}(5), 1287-1294.
#' \doi{10.2307/1911963}
#'
#' Goldfeld SM, Quandt RE (1965). "Some tests for homoscedasticity."
#' \emph{Journal of the American Statistical Association}, \bold{60}(310),
#' 539-547. \doi{10.1080/01621459.1965.10480811}
#'
#' Koenker R (1981). "A note on studentizing a test for heteroscedasticity."
#' \emph{Journal of Econometrics}, \bold{17}(1), 107-112.
#' \doi{10.1016/0304-4076(81)90062-2}
#'
#' Li Z, Yao J (2019). "Testing for heteroscedasticity in high-dimensional
#' regressions." \emph{Econometrics and Statistics}, \bold{9}, 122-139.
#' \doi{10.1016/j.ecosta.2018.01.001}
#'
#' Rana MS, Midi H, Imon AHMR (2008). "A robust modification of the
#' Goldfeld-Quandt test for the detection of heteroscedasticity in the presence
#' of outliers." \emph{Journal of Mathematics and Statistics}, \bold{4}(4),
#' 277-283. \doi{10.3844/jmssp.2008.277.283}
#'
#' White H (1980). "A heteroskedasticity-consistent covariance matrix estimator
#' and a direct test for heteroskedasticity." \emph{Econometrica},
#' \bold{48}(4), 817-838. \doi{10.2307/1912934}
#'
#' Wilcox RR, Keselman HJ (2006). "Detecting heteroscedasticity in a simple
#' regression model via quantile regression slopes." \emph{Journal of
#' Statistical Computation and Simulation}, \bold{76}(8), 705-712.
#' \doi{10.1080/10629360500107923}
#'
#' Zhou QM, Song PX-K, Thompson ME (2015). "Profiling heteroscedasticity in
#' linear regression models." \emph{Canadian Journal of Statistics},
#' \bold{43}(3), 358-377. \doi{10.1002/cjs.11252}
#'
#' @seealso \code{\link{kah3.test}}, \code{\link{kah.robust.test}}.
#'
#' @examples
#' run.all.het(dist ~ speed, data = cars, include.optional = FALSE)
#'
#' @export
run.all.het <- function(x, data = NULL, order.by = NULL, alpha = 0.75,
                        include.optional = TRUE, bootstrap = FALSE, B = 499,
                        seed = 1) {
  d <- .het_prepare(x, data, order.by)
  y <- d$y; Z <- d$Z; ord <- d$ord
  fit <- stats::lm.fit(cbind(1, Z), y)
  e <- fit$residuals

  rows <- list()
  add <- function(test, group, expr) {
    r <- tryCatch(expr, error = function(err) list(stat = NA_real_, df = NA_character_,
                                                     p = NA_real_, note = conditionMessage(err)))
    rows[[length(rows) + 1]] <<- data.frame(Test = test, Group = group,
                                            Statistic = unname(r$stat), df = as.character(r$df),
                                            p_value = unname(r$p),
                                            Note = if (is.null(r$note)) "" else r$note,
                                            stringsAsFactors = FALSE)
  }
  from_htest <- function(h, df) list(stat = h$statistic, df = df, p = h$p.value)

  add("KaH-III", "Proposed", {
    h <- .kah3_core(d)
    from_htest(h, sprintf("3, %g", h$parameter["df"]))
  })
  add(sprintf("KaH robust (alpha = %s)", alpha), "Proposed", {
    h <- .kah_robust_core(d, alpha, "fmax", B, seed)
    from_htest(h, sprintf("3, %.2f", h$parameter["df"]))
  })
  if (bootstrap) add(sprintf("KaH robust bootstrap (alpha = %s)", alpha), "Proposed", {
    h <- .kah_robust_core(d, alpha, "bootstrap", B, seed)
    from_htest(h, sprintf("B = %d", B))
  })
  add("Goldfeld-Quandt", "Classical", .gq_test(Z, y, ord))
  add("Breusch-Pagan (Koenker)", "Classical", .bp_koenker(Z, e, Z))
  add("White", "Classical", .bp_koenker(Z, e, .white_terms(Z)))
  add("MGQ (Rana et al. 2008)", "Robust", .with_seed(seed, .mgq_test(Z, y, ord)))

  if (include.optional) {
    if (requireNamespace("skedastic", quietly = TRUE)) {
      Zd <- as.data.frame(Z); names(Zd) <- paste0("x", seq_len(ncol(Z)))   # safe names
      m0 <- stats::lm(y ~ ., data = cbind(y = y, Zd))
      sk <- function(fun, ...) {
        h <- suppressWarnings(fun(m0, ...))
        list(stat = h$statistic, df = "", p = h$p.value)
      }
      add("Wilcox-Keselman (2006)", "Robust", .with_seed(seed, sk(skedastic::wilcox_keselman)))
      add("Zhou-Song-Thompson (2015)", "Recent", sk(skedastic::zhou_etal, method = "pooled", seed = seed))
      add("Li-Yao (2019)", "Recent", sk(skedastic::li_yao, method = "cvt"))
    } else {
      message("Install 'skedastic' to add the Wilcox-Keselman, Zhou-Song-Thompson and Li-Yao tests.")
    }
  }

  out <- do.call(rbind, rows)
  attr(out, "data.name") <- d$data.name
  attr(out, "order.by") <- d$ord.name
  class(out) <- c("het_battery", "data.frame")
  out
}

#' @export
print.het_battery <- function(x, digits = 4, ...) {
  cat("\nHeteroscedasticity tests for:", attr(x, "data.name"), "\n")
  cat("Data sorted by:", attr(x, "order.by"), "\n\n")
  df <- as.data.frame(unclass(x))[, c("Test", "Group", "Statistic", "df", "p_value")]
  df$Statistic <- formatC(df$Statistic, digits = digits, format = "g")
  df$p_value <- ifelse(is.na(x$p_value), "NA", format.pval(x$p_value, digits = digits, eps = 1e-4))
  df$` ` <- ifelse(!is.na(x$p_value) & x$p_value < 0.05, "*", "")
  for (g in unique(df$Group)) {
    cat(g, "\n")
    print(df[df$Group == g, setdiff(names(df), "Group")], row.names = FALSE, right = FALSE)
    cat("\n")
  }
  notes <- x$Note[nzchar(x$Note)]
  if (length(notes)) cat("Notes:", paste0(x$Test[nzchar(x$Note)], ": ", notes, collapse = "; "), "\n")
  cat("* p < 0.05\n")
  invisible(x)
}

# ---- in-house implementations of the comparison tests ----------------------

.gq_test <- function(Z, y, ord) {
  o <- order(ord); y <- y[o]; Z <- Z[o, , drop = FALSE]
  n <- length(y); p <- ncol(Z); m <- floor(n / 3)
  g1 <- 1:m; g2 <- (n - m + 1):n
  rss <- function(i) sum(stats::lm.fit(cbind(1, Z[i, , drop = FALSE]), y[i])$residuals^2)
  df <- m - p - 1
  f <- (rss(g2) / df) / (rss(g1) / df)
  list(stat = f, df = sprintf("%d, %d", df, df),
       p = 2 * min(stats::pf(f, df, df), stats::pf(f, df, df, lower.tail = FALSE)))
}

.bp_koenker <- function(Z, e, W) {
  n <- length(e); u <- e^2
  r2 <- 1 - sum(stats::lm.fit(cbind(1, W), u)$residuals^2) / sum((u - mean(u))^2)
  stat <- n * r2
  list(stat = stat, df = ncol(W), p = stats::pchisq(stat, ncol(W), lower.tail = FALSE))
}

.white_terms <- function(Z) {
  W <- cbind(Z, Z^2)
  if (ncol(Z) > 1) for (i in 1:(ncol(Z) - 1)) for (j in (i + 1):ncol(Z)) W <- cbind(W, Z[, i] * Z[, j])
  W
}

.mgq_test <- function(Z, y, ord) {
  o <- order(ord); y <- y[o]; Z <- Z[o, , drop = FALSE]
  n <- length(y); p <- ncol(Z); m <- floor(n / 3)
  lts <- robustbase::ltsReg(x = Z, y = y, intercept = TRUE, mcd = FALSE)
  clean <- lts$lts.wt == 1
  X <- cbind(1, Z)
  b <- stats::lm.fit(X[clean, , drop = FALSE], y[clean])$coefficients
  r <- drop(y - X %*% b)
  Xc <- X[clean, , drop = FALSE]
  h <- rowSums((Xc %*% solve(crossprod(Xc))) * Xc)
  r[clean] <- r[clean] / (1 - h)
  stat <- stats::median(r[(n - m + 1):n]^2) / stats::median(r[1:m]^2)
  df <- m - p - 1
  list(stat = stat, df = sprintf("%d, %d", df, df),
       p = 2 * min(stats::pf(stat, df, df), stats::pf(stat, df, df, lower.tail = FALSE)))
}
