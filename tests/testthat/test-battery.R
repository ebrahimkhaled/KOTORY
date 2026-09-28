test_that("run.all.het returns one row per test and never aborts", {
  b <- run.all.het(dist ~ speed, data = cars, include.optional = FALSE)
  expect_s3_class(b, "het_battery")
  expect_equal(nrow(b), 6)
  expect_true(all(b$p_value >= 0 & b$p_value <= 1))
  expect_output(print(b), "Proposed")
})

test_that("the in-house classical tests agree with lmtest", {
  skip_if_not_installed("lmtest")
  fit <- lm(dist ~ speed, data = cars)
  b <- run.all.het(fit, include.optional = FALSE)
  expect_equal(b$p_value[b$Test == "Breusch-Pagan (Koenker)"], unname(lmtest::bptest(fit)$p.value))
  w <- lmtest::bptest(fit, ~ speed + I(speed^2), data = cars)
  expect_equal(b$p_value[b$Test == "White"], unname(w$p.value))
  gq <- lmtest::gqtest(fit, order.by = ~ speed, data = cars, fraction = 50 - 2 * 16,
                       alternative = "two.sided")
  expect_equal(b$p_value[b$Test == "Goldfeld-Quandt"], unname(gq$p.value))
})

test_that("non-syntactic regressor names work", {
  b <- run.all.het(dist ~ log(speed), data = cars, include.optional = FALSE)
  expect_false(anyNA(b$p_value))
})

test_that("optional tests run when skedastic is installed", {
  skip_on_cran()
  skip_if_not_installed("skedastic")
  b <- run.all.het(dist ~ speed, data = cars)
  expect_equal(nrow(b), 9)
})
