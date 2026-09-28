# Internal helpers shared by the tests. Not exported.

# Accept an lm object or a formula (+ data); return the response, the
# regressor matrix without the intercept column, and the sorting variable.
.het_prepare <- function(x, data = NULL, order.by = NULL, caller = "test") {
  if (inherits(x, "lm")) {
    mf <- stats::model.frame(x)
    X <- stats::model.matrix(x)
    y <- stats::model.response(mf)
    dname <- paste(deparse(stats::formula(x)), collapse = " ")
  } else if (inherits(x, "formula")) {
    if (is.null(data)) data <- environment(x)
    mf <- stats::model.frame(x, data = data, na.action = stats::na.omit)
    X <- stats::model.matrix(x, mf)
    y <- stats::model.response(mf)
    dname <- paste(deparse(x), collapse = " ")
  } else {
    stop("'x' must be a fitted 'lm' object or a model formula.", call. = FALSE)
  }
  if (!is.numeric(y)) stop("The response must be numeric.", call. = FALSE)
  int <- colnames(X) == "(Intercept)"
  if (!any(int)) stop("The model must contain an intercept.", call. = FALSE)
  Z <- X[, !int, drop = FALSE]
  if (ncol(Z) == 0) stop("The model must contain at least one regressor.", call. = FALSE)

  ord <- .het_order(order.by, Z, mf, data, nrow(Z))
  list(y = as.numeric(y), Z = Z, ord = ord$values, ord.name = ord$name,
       n = length(y), p = ncol(Z), data.name = dname)
}

# The sorting variable: default is the first regressor (as in the thesis);
# a regressor name, a column of 'data', or a numeric vector are accepted.
.het_order <- function(order.by, Z, mf, data, n) {
  if (is.null(order.by)) return(list(values = Z[, 1], name = colnames(Z)[1]))
  if (is.character(order.by) && length(order.by) == 1) {
    if (order.by %in% colnames(Z)) return(list(values = Z[, order.by], name = order.by))
    if (order.by %in% names(mf)) return(list(values = mf[[order.by]], name = order.by))
    if (!is.null(data) && is.data.frame(data) && order.by %in% names(data) && nrow(data) == n)
      return(list(values = data[[order.by]], name = order.by))
    stop("'order.by' = \"", order.by, "\" was not found among the regressors or in 'data'.", call. = FALSE)
  }
  if (is.numeric(order.by) && length(order.by) == n)
    return(list(values = order.by, name = "user-supplied variable"))
  stop("'order.by' must be a regressor name or a numeric vector of length n.", call. = FALSE)
}

# The three parts exactly as in the thesis code: each of size floor(n/3);
# when n is not a multiple of 3 the one or two observations at the part
# boundaries are left out.
.kah_parts <- function(n) {
  list(1:floor(n / 3),
       (round(n / 3) + 1):floor(2 * n / 3),
       ceiling(2 * n / 3 + 1):n)
}

# Run code under a fixed seed without disturbing the caller's RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  expr
}

# Raw LTS scale and standardized raw residuals of one part. Small parts are searched
# exhaustively (all elemental subsets), which is deterministic and also succeeds when
# repeated regressor values leave few non-singular random subsets.
.lts_part <- function(Z, y, alpha) {
  m <- length(y)
  nsamp <- if (choose(m, ncol(Z) + 1) <= 3000) "exact" else 500
  fit <- tryCatch(
    robustbase::ltsReg(x = Z, y = y, intercept = TRUE, alpha = alpha, mcd = FALSE, nsamp = nsamp),
    error = function(e) stop(sprintf(paste0(
      "The LTS fit failed in a third of the data with %d observations (%s). ",
      "This usually means the regressors take too few distinct values within that third."),
      m, conditionMessage(e)), call. = FALSE))
  r <- drop(y - cbind(1, Z) %*% fit$raw.coefficients) / fit$raw.scale
  list(scale = fit$raw.scale, resid = r)
}

.lts3 <- function(Z, y, alpha) {
  parts <- .kah_parts(length(y))
  res <- lapply(parts, function(i) .lts_part(Z[i, , drop = FALSE], y[i], alpha))
  list(scale = vapply(res, `[[`, numeric(1), "scale"),
       resid = unlist(lapply(res, `[[`, "resid")))
}

.check_n <- function(m, p, robust) {
  need <- if (robust) 2 * (p + 1) + 1 else p + 2
  if (m < need)
    stop(sprintf("Each third has %d observations; at least %d are needed for %s with %d regressor(s).",
                 m, need, if (robust) "the LTS fit" else "the OLS fit", p), call. = FALSE)
}
