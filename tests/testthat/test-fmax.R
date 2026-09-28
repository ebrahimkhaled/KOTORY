test_that("qfmax reproduces Hartley's published values (k = 3)", {
  # Pearson & Hartley, Table 31: df = 4 -> 15.5 (5%), 37 (1%); df = 10 -> 4.85, 7.4
  expect_equal(qfmax(0.95, 4), 15.5, tolerance = 0.01)
  expect_equal(qfmax(0.99, 4), 37, tolerance = 0.02)
  expect_equal(qfmax(0.95, 10), 4.85, tolerance = 0.01)
  expect_equal(qfmax(0.99, 10), 7.4, tolerance = 0.02)
})

test_that("k = 2 equals the two-sided F distribution", {
  for (v in c(3, 8, 25)) for (cc in c(1.5, 3, 7)) {
    expect_equal(pfmax(cc, v, k = 2), pf(cc, v, v) - pf(1 / cc, v, v), tolerance = 1e-7)
  }
})

test_that("pfmax and qfmax are inverse and the tails add to one", {
  p <- c(0.5, 0.9, 0.95, 0.99, 0.995)
  expect_equal(pfmax(qfmax(p, 7.3), 7.3), p, tolerance = 1e-7)
  expect_equal(pfmax(3, 6) + pfmax(3, 6, lower.tail = FALSE), 1)
  expect_equal(qfmax(0.05, 6, lower.tail = FALSE), qfmax(0.95, 6))
})

test_that("dfmax integrates to one and is the derivative of pfmax", {
  expect_equal(integrate(function(x) dfmax(x, 6), 1, Inf)$value, 1, tolerance = 1e-5)
  h <- 1e-5
  expect_equal(dfmax(2.5, 9), (pfmax(2.5 + h, 9) - pfmax(2.5 - h, 9)) / (2 * h), tolerance = 1e-5)
})

test_that("rfmax agrees with pfmax", {
  set.seed(1)
  r <- rfmax(4000, 5)
  expect_gt(suppressWarnings(ks.test(r, function(q) pfmax(q, 5))$p.value), 0.001)
})

test_that("edge cases and argument checks", {
  expect_equal(pfmax(c(0.5, 1), 5), c(0, 0))
  expect_equal(pfmax(Inf, 5), 1)
  expect_equal(qfmax(c(0, 1), 5), c(1, Inf))
  expect_true(is.na(pfmax(NA, 5)))
  expect_error(pfmax(2, -1), "positive")
  expect_error(pfmax(2, 5, k = 1.5), "integer")
})
