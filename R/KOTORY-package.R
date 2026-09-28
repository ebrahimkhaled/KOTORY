#' KOTORY: Robust Three-Group Tests for Heteroscedasticity in Linear Regression
#'
#' Tests for heteroscedasticity that sort the data by a regressor, split them
#' into three equal parts and compare the error scale of the parts.
#'
#' @section Choosing a test:
#' \describe{
#'   \item{Clean data, normal errors}{\code{\link{kah3.test}} --- least squares
#'     in each part; its null distribution is Hartley's maximum F-ratio with
#'     three groups, exactly.}
#'   \item{Outliers possible}{\code{\link{kah.robust.test}} --- least trimmed
#'     squares in each part, so outliers neither create nor hide
#'     heteroscedasticity as easily.}
#'   \item{Errors possibly non-normal}{\code{kah.robust.test(method =
#'     "bootstrap")} --- the reference is resampled from the data's own
#'     residuals.}
#'   \item{Comparing with other tests}{\code{\link{run.all.het}} --- the KaH
#'     tests next to Goldfeld-Quandt, Breusch-Pagan, White, the robust MGQ and,
#'     when \pkg{skedastic} is installed, three more recent tests.}
#' }
#'
#' @section Distribution functions:
#' \code{\link{pfmax}}, \code{\link{dfmax}}, \code{\link{qfmax}} and
#' \code{\link{rfmax}} give Hartley's maximum F-ratio for any number of groups
#' and any positive degrees of freedom.
#'
#' @section Origin:
#' The three-group tests were proposed in the doctoral thesis of Ahmed
#' El-Kotory (Alexandria University), where their critical values were
#' tabulated by simulation. This package replaces the tables by the exact
#' distribution (least squares version) and by an effective-degrees-of-freedom
#' approximation, a Monte Carlo or a bootstrap reference (robust version).
#'
#' @keywords internal
#' @concept heteroscedasticity
#' @concept linear regression
#' @concept robust statistics
#' @concept least trimmed squares
#' @concept Goldfeld-Quandt
#' @concept Hartley Fmax
"_PACKAGE"
