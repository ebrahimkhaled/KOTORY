sim_data <- function(n = 60, p = 1, seed = 3, het = 0) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  y <- drop(1 + X %*% rep(1, p) + rnorm(n) * exp(het * X[, 1]))
  data.frame(y = y, X)
}

test_that("kah3.test matches a hand computation of the thesis statistic", {
  d <- sim_data(62, 2)
  h <- kah3.test(y ~ x1 + x2, data = d)
  s <- d[order(d$x1), ]
  parts <- list(1:20, 22:41, 43:62)          # thesis split for n = 62
  ms <- sapply(parts, function(i) anova(lm(y ~ x1 + x2, data = s[i, ]))[3, 3])
  expect_equal(unname(h$statistic), max(ms) / min(ms))
  expect_equal(unname(h$parameter["df"]), 20 - 3)
  expect_equal(h$p.value, pfmax(max(ms) / min(ms), 17, lower.tail = FALSE))
  expect_equal(h$n.used, 60)
  expect_s3_class(h, "htest")
})

test_that("formula and lm input agree; statistic is scale invariant", {
  d <- sim_data()
  a <- kah3.test(y ~ x1, data = d)
  b <- kah3.test(lm(y ~ x1, data = d))
  expect_equal(a$statistic, b$statistic)
  d2 <- d; d2$y <- 10 * d2$y + 3
  expect_equal(unname(kah3.test(y ~ x1, data = d2)$statistic), unname(a$statistic))
})

test_that("order.by accepts a name or a vector", {
  d <- sim_data(60, 2)
  a <- kah3.test(y ~ x1 + x2, data = d, order.by = "x2")
  b <- kah3.test(y ~ x1 + x2, data = d, order.by = d$x2)
  expect_equal(a$statistic, b$statistic)
  expect_error(kah3.test(y ~ x1, data = d, order.by = "nope"), "not found")
})

test_that("the critical values are the Fmax quantiles", {
  h <- kah3.test(y ~ x1, data = sim_data())
  expect_equal(unname(h$critical["5%"]), qfmax(0.95, h$parameter["df"]))
})

test_that("kah.robust.test matches a direct LTS computation", {
  d <- sim_data(60)
  h <- kah.robust.test(y ~ x1, data = d, alpha = 0.9, seed = 1)
  s <- d[order(d$x1), ]
  sc <- KOTORY:::.with_seed(1, sapply(list(1:20, 21:40, 41:60), function(i)
    robustbase::ltsReg(x = as.matrix(s$x1[i]), y = s$y[i], alpha = 0.9, mcd = FALSE)$raw.scale))
  expect_equal(unname(h$statistic), max(sc) / min(sc))
  expect_equal(h$p.value, pfmax(max(sc)^2 / min(sc)^2, kah.nu.star(20, 1, 0.9), lower.tail = FALSE))
})

test_that("the robust test does not disturb the user's random numbers", {
  d <- sim_data()
  set.seed(99); a <- runif(1)
  set.seed(99); invisible(kah.robust.test(y ~ x1, data = d)); b <- runif(1)
  expect_equal(a, b)
})

test_that("mc and bootstrap references run and are reproducible", {
  d <- sim_data(45)
  a <- kah.robust.test(y ~ x1, data = d, method = "mc", B = 39)
  b <- kah.robust.test(y ~ x1, data = d, method = "mc", B = 39)
  expect_equal(a$p.value, b$p.value)
  bs <- kah.robust.test(y ~ x1, data = d, method = "bootstrap", B = 39)
  expect_true(bs$p.value > 0 && bs$p.value <= 1)
  expect_false("df" %in% names(bs$parameter))
})

test_that("an outlier fools the least squares test but not the robust one", {
  set.seed(2)
  x <- runif(60); y0 <- 1 + x + rnorm(60); y <- y0; y[which.min(x)] <- 12
  expect_gt(kah3.test(y0 ~ x)$p.value, 0.1)
  expect_lt(kah3.test(y ~ x)$p.value, 0.001)
  # the LTS fits trim the outlier: the robust statistic does not move at all
  expect_equal(kah.robust.test(y ~ x)$statistic, kah.robust.test(y0 ~ x)$statistic)
})

test_that("strong heteroscedasticity is detected", {
  d <- sim_data(150, het = 2)
  expect_lt(kah3.test(y ~ x1, data = d)$p.value, 0.01)
  expect_lt(kah.robust.test(y ~ x1, data = d)$p.value, 0.01)
})

test_that("repeated regressor values (grouped design) work", {
  # Pindyck-Rubinfeld housing data: x takes 4 values, 5 times each
  x <- rep(c(5, 10, 15, 20), each = 5)
  y <- c(1.8, 2, 2, 2, 2.1, 3.1, 3.2, 3.5, 3.5, 3.6, 4.2, 4.2, 4.5, 4.8, 5, 4.8, 5, 5.7, 6, 6.2)
  h <- kah.robust.test(y ~ x)
  expect_true(is.finite(h$statistic) && h$p.value >= 0 && h$p.value <= 1)
})

test_that("input checks", {
  d <- sim_data(12, 3)
  expect_error(kah.robust.test(y ~ x1 + x2 + x3, data = d), "at least")
  expect_error(kah3.test(y ~ x1 - 1, data = sim_data()), "intercept")
  expect_error(kah.robust.test(y ~ x1, data = sim_data(), alpha = 0.3), "between")
  expect_error(kah.nu.star(20, 1, 0.6), "mc")
  expect_error(kah3.test("y ~ x"), "formula")
})
