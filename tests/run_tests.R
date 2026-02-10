source("postSim.R")
source("post.R")
library(testthat)

# Existing test: post rejects silently recycled glm weights
test_that("post rejects silently recycled glm weights", {
  fit <- glm(vs ~ wt + hp, data = mtcars, family = binomial())
  expect_error(
    post(fit, weights = rep(1, 2), n.sims = 20),
    "weights must have the same length as the estimation sample"
  )
})

# New test: k is not clobbered by holds loop in glm x1+x2 branch
test_that("k variable not shadowed by holds in glm x1+x2 branch", {
  fit <- glm(vs ~ wt + hp + am, data = mtcars, family = binomial())
  set.seed(42)
  res <- post(fit,
              x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(am = 1),
              n.sims = 50)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

# New test: multiple holds don't corrupt k
test_that("multiple holds entries don't corrupt k in x1+x2 branch", {
  fit <- glm(vs ~ wt + hp + am + gear, data = mtcars, family = binomial())
  set.seed(42)
  res <- post(fit,
              x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(am = 1, gear = 4),
              n.sims = 50)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

cat("\nAll tests passed!\n")
