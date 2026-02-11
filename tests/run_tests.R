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

# Factor x1vals: factor defined outside formula (as a column in data)
test_that("factor x1vals works when factor defined outside formula", {
  d <- mtcars
  d$am_f <- factor(ifelse(d$am == 1, "manual", "auto"))
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  set.seed(1)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = 50)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# Factor x1vals: factor defined inside the model call via transform()
test_that("factor x1vals works when factor defined inside model call", {
  fit <- glm(vs ~ am_f + wt,
             data = transform(mtcars, am_f = factor(ifelse(am == 1, "manual", "auto"))),
             family = binomial())
  set.seed(1)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = 50)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# Invalid factor level produces informative error
test_that("invalid factor x1vals produces error listing valid levels", {
  d <- mtcars
  d$am_f <- factor(ifelse(d$am == 1, "manual", "auto"))
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "am_f", x1vals = c("auto", "INVALID"), n.sims = 20),
    "x1vals must be valid levels of factor 'am_f'"
  )
})

# Numeric-like factor levels work when passed as numbers
test_that("numeric-like factor levels are matched after character conversion", {
  d <- mtcars
  d$cyl_f <- factor(d$cyl)
  fit <- glm(vs ~ cyl_f + wt, data = d, family = binomial())
  set.seed(1)
  res <- post(fit, x1name = "cyl_f", x1vals = c(4, 6, 8), n.sims = 50)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# Inline formula transformation produces informative error
test_that("inline formula transformation in x1name gives clear error", {
  fit <- glm(vs ~ factor(am) + wt, data = mtcars, family = binomial())
  expect_error(
    post(fit, x1name = "factor(am)", x1vals = c("0", "1"), n.sims = 20),
    "appears to use an inline transformation"
  )
})

# Misspelled x1name produces informative error
test_that("misspelled x1name produces error naming the variable", {
  fit <- glm(vs ~ wt + hp, data = mtcars, family = binomial())
  expect_error(
    post(fit, x1name = "weight", x1vals = c(2, 4), n.sims = 20),
    "x1name='weight' is not a variable in the model"
  )
})

# Misspelled x2name produces informative error
test_that("misspelled x2name produces error naming the variable", {
  fit <- glm(vs ~ wt + hp, data = mtcars, family = binomial())
  expect_error(
    post(fit, x1name = "wt", x1vals = c(2, 4),
         x2name = "horsepower", x2vals = c(100, 200), n.sims = 20),
    "x2name='horsepower' is not a variable in the model"
  )
})

cat("\nAll tests passed!\n")
