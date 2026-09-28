# Builds the internal table of effective degrees of freedom (R/sysdata.rda)
# used by kah.nu.star(). One part of m observations, p regressors, normal
# errors; nu* solves trigamma(nu*/2) = var(log raw LTS scale^2).
# Each grid cell is cached in data-raw/nu_cells/ so the run can resume.

library(parallel)
cells_dir <- "data-raw/nu_cells"; dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)
R <- 20000
grid <- expand.grid(m = c(5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50, 70, 100, 150, 200),
                    p = 1:5, alpha = c(0.5, 0.75, 0.9))
grid <- grid[grid$m >= 2 * (grid$p + 1) + 1, ]

cell <- function(k) {
  m <- grid$m[k]; p <- grid$p[k]; a <- grid$alpha[k]
  f <- file.path(cells_dir, sprintf("a%.2f_p%d_m%03d.rds", a, p, m))
  if (file.exists(f)) return(readRDS(f))
  set.seed(1000 * m + 10 * p + round(100 * a))
  s2 <- vapply(seq_len(R), function(r) {
    Z <- scale(matrix(stats::runif(3 * m * p), 3 * m, p))[1:m, , drop = FALSE]
    Z <- Z[order(Z[, 1]), , drop = FALSE]
    y <- stats::rnorm(m)
    robustbase::ltsReg(x = Z, y = y, intercept = TRUE, alpha = a, mcd = FALSE)$raw.scale^2
  }, numeric(1))
  v <- stats::var(log(s2[s2 > 0]))
  out <- data.frame(alpha = a, p = p, m = m, nu = m - p - 1,
                    nu_star = 2 * stats::uniroot(function(z) trigamma(z) - v, c(1e-3, 1e5))$root,
                    zero_frac = mean(s2 <= 0))
  saveRDS(out, f)
  out
}

cl <- makeCluster(max(1, detectCores() - 2))
clusterExport(cl, c("grid", "R", "cells_dir", "cell"))
res <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(grid)), function(k) tryCatch(cell(k), error = function(e) NULL)))
stopCluster(cl)

res$ratio <- res$nu_star / res$nu
.nu_table <- res[order(res$alpha, res$p, res$m), c("alpha", "p", "m", "nu", "nu_star", "ratio")]
rownames(.nu_table) <- NULL
print(.nu_table, digits = 3)
save(.nu_table, file = "R/sysdata.rda", compress = "xz")
