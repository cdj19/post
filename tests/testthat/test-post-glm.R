if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required")
}

testthat::test_that("holds preserves non-numeric types for glm models", {
  source("postSim.R", local = TRUE)
  source("post.R", local = TRUE)

  d <- mtcars
  d$am_f <- factor(ifelse(d$am == 1, "manual", "auto"))
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())

  set.seed(1)
  res <- post(fit, holds = list(am_f = "manual"), n.sims = 50, digits = 6)

  testthat::expect_s4_class(res, "post")
  testthat::expect_false(any(is.na(res@est)))
})

testthat::test_that("post rejects silently recycled glm weights", {
  source("postSim.R", local = TRUE)
  source("post.R", local = TRUE)

  fit <- glm(vs ~ wt + hp, data = mtcars, family = binomial())

  testthat::expect_error(
    post(fit, weights = rep(1, 2), n.sims = 20),
    "weights must have the same length as the estimation sample"
  )
})
