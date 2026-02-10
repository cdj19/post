if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required")
}

testthat::test_that("post rejects silently recycled polr weights", {
  testthat::skip_if_not_installed("MASS")

  source("postSim.R", local = TRUE)
  source("post.R", local = TRUE)

  d <- transform(mtcars, gear_f = ordered(gear))
  fit <- MASS::polr(gear_f ~ wt + hp, data = d, method = "logistic")

  testthat::expect_error(
    post(fit, weights = rep(1, 3), n.sims = 20),
    "weights must have the same length as the estimation sample"
  )
})
