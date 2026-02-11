source("postSim.R")
source("post.R")
library(testthat)
library(MASS)

set.seed(12345)

# ── Data setup ──────────────────────────────────────────────────────────────

d <- mtcars
d$am_f  <- factor(ifelse(d$am == 1, "manual", "auto"))
d$vs_f  <- factor(ifelse(d$vs == 1, "straight", "vee"))
d$cyl_f <- factor(d$cyl)
d$gear_o <- ordered(d$gear)   # levels: "3","4","5"
d$cyl_o  <- ordered(d$cyl)    # levels: "4","6","8"

N <- 100  # n.sims – fast but enough to exercise all code paths

# ============================================================================
# 1. NO FACTOR VARIABLES, NO INTERACTIONS
# ============================================================================

# ── lm ─────────────────────────────────────────────────────────────────────
test_that("lm: no x (marginal prediction)", {
  fit <- lm(mpg ~ wt + hp, data = d)
  res <- post(fit, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: x1 only", {
  fit <- lm(mpg ~ wt + hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_equal(nrow(res@est), 3)  # 2 values + first difference
})

test_that("lm: x1 + x2", {
  fit <- lm(mpg ~ wt + hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("lm: with holds", {
  fit <- lm(mpg ~ wt + hp + drat, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              holds = list(drat = 3.5), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm (logit) ───────────────────────────────────────────────────────────
test_that("glm logit: no x", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_true(res@est[1, 1] >= 0 && res@est[1, 1] <= 1)
})

test_that("glm logit: x1 only", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 3, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_equal(nrow(res@est), 4)  # 3 vals + first diff
})

test_that("glm logit: x1 + x2", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("glm logit: x1 + x2 with holds", {
  fit <- glm(vs ~ wt + hp + drat, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(drat = 3.5), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm (probit) ──────────────────────────────────────────────────────────
test_that("glm probit: x1 only", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial(link = "probit"))
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm (cloglog) ─────────────────────────────────────────────────────────
test_that("glm cloglog: x1 only", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial(link = "cloglog"))
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm (poisson / log link) ──────────────────────────────────────────────
test_that("glm poisson: x1 only", {
  fit <- glm(carb ~ wt + hp, data = d, family = poisson())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm poisson: x1 + x2", {
  fit <- glm(carb ~ wt + hp, data = d, family = poisson())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr (logistic) ───────────────────────────────────────────────────────
test_that("polr logistic: no x", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_equal(nrow(res@est), length(levels(d$gear_o)))
})

test_that("polr logistic: x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr logistic: x1 + x2", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("polr logistic: cut, no x", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_true(res@est[1, 1] >= 0 && res@est[1, 1] <= 1)
})

test_that("polr logistic: cut, x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr logistic: cut, x1 + x2", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

# ── polr (probit) ─────────────────────────────────────────────────────────
test_that("polr probit: x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "probit")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr probit: cut, x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "probit")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr (cloglog) ────────────────────────────────────────────────────────
test_that("polr cloglog: no x", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: x1 + x2", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: cut, no x", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_true(res@est[1, 1] >= 0 && res@est[1, 1] <= 1)
})

test_that("polr cloglog: cut, x1 only", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: cut, x1 + x2", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})


# ============================================================================
# 2. FACTOR VARIABLES, NO INTERACTIONS
# ============================================================================

# ── lm ─────────────────────────────────────────────────────────────────────
test_that("lm: factor x1, no interactions", {
  fit <- lm(mpg ~ cyl_f + wt, data = d)
  res <- post(fit, x1name = "cyl_f", x1vals = c(4, 6, 8), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: factor x1 + factor x2, no interactions", {
  fit <- lm(mpg ~ am_f + vs_f + wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: factor x1 + numeric x2, no interactions", {
  fit <- lm(mpg ~ am_f + wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: factor in holds", {
  fit <- lm(mpg ~ am_f + wt + hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              holds = list(am_f = "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm ────────────────────────────────────────────────────────────────────
test_that("glm: factor x1, no interactions", {
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm: 3-level factor x1, no interactions", {
  fit <- glm(vs ~ cyl_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "cyl_f", x1vals = c(4, 6, 8), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_equal(nrow(res@est), 4)  # 3 levels + first diff
})

test_that("glm: factor x1 + numeric x2, no interactions", {
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr ───────────────────────────────────────────────────────────────────
test_that("polr: factor x1, no interactions", {
  fit <- polr(gear_o ~ am_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor x1, cut", {
  fit <- polr(gear_o ~ am_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor x1 + numeric x2, no interactions", {
  fit <- polr(gear_o ~ am_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor x1 + numeric x2, cut", {
  fit <- polr(gear_o ~ am_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor in holds", {
  fit <- polr(gear_o ~ am_f + wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              holds = list(am_f = "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})


# ============================================================================
# 3. INTERACTIONS, NO FACTOR VARIABLES
# ============================================================================

# ── lm ─────────────────────────────────────────────────────────────────────
test_that("lm: numeric interaction, x1 only", {
  fit <- lm(mpg ~ wt * hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: numeric interaction, x1 + x2 (the interacted vars)", {
  fit <- lm(mpg ~ wt * hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("lm: numeric interaction, x1 + x2 + holds", {
  fit <- lm(mpg ~ wt * hp + drat, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(drat = 3.5), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm ────────────────────────────────────────────────────────────────────
test_that("glm: numeric interaction, x1 only", {
  fit <- glm(vs ~ wt * hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm: numeric interaction, x1 + x2", {
  fit <- glm(vs ~ wt * hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr ───────────────────────────────────────────────────────────────────
test_that("polr: numeric interaction, x1 only", {
  fit <- polr(gear_o ~ wt * hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: numeric interaction, x1 + x2", {
  fit <- polr(gear_o ~ wt * hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: numeric interaction, cut, x1 + x2", {
  fit <- polr(gear_o ~ wt * hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: numeric interaction, cut, x1 + x2", {
  fit <- polr(gear_o ~ wt * am, data = d, method = "cloglog")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "am", x2vals = c(0, 1), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})


# ============================================================================
# 4. FACTOR × NUMERIC INTERACTIONS
# ============================================================================

# ── lm ─────────────────────────────────────────────────────────────────────
test_that("lm: factor*numeric, factor as x1", {
  fit <- lm(mpg ~ am_f * wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: factor*numeric, factor x1 + numeric x2", {
  fit <- lm(mpg ~ am_f * wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("lm: factor*numeric, numeric x1 + factor x2", {
  fit <- lm(mpg ~ am_f * wt, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "am_f", x2vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: 3-level factor * numeric, factor x1 + numeric x2", {
  fit <- lm(mpg ~ cyl_f * wt, data = d)
  res <- post(fit, x1name = "cyl_f", x1vals = c("4", "6", "8"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm ────────────────────────────────────────────────────────────────────
test_that("glm: factor*numeric, factor x1 + numeric x2", {
  fit <- glm(vs ~ am_f * wt, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm: factor*numeric, numeric x1 + factor x2", {
  fit <- glm(vs ~ am_f * wt, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "am_f", x2vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm: factor*numeric with holds on a third variable", {
  fit <- glm(vs ~ am_f * wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4),
              holds = list(hp = 150), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr ───────────────────────────────────────────────────────────────────
test_that("polr: factor*numeric, factor x1 + numeric x2", {
  fit <- polr(gear_o ~ am_f * wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor*numeric, numeric x1 + factor x2", {
  fit <- polr(gear_o ~ am_f * wt, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "am_f", x2vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor*numeric, cut, factor x1 + numeric x2", {
  fit <- polr(gear_o ~ am_f * wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: factor*numeric, cut, x1 + x2", {
  fit <- polr(gear_o ~ am_f * wt, data = d, method = "cloglog")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})


# ============================================================================
# 5. FACTOR × FACTOR INTERACTIONS
# ============================================================================

# ── lm ─────────────────────────────────────────────────────────────────────
test_that("lm: factor*factor, x1 + x2", {
  fit <- lm(mpg ~ am_f * vs_f + wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("lm: factor*factor, x1 only (interaction still in model)", {
  fit <- lm(mpg ~ am_f * vs_f + wt, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: 3-level * 2-level factor interaction, x1 + x2", {
  fit <- lm(mpg ~ cyl_f * am_f + wt, data = d)
  res <- post(fit, x1name = "cyl_f", x1vals = c("4", "6", "8"),
              x2name = "am_f", x2vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("lm: factor*factor with holds on third variable", {
  fit <- lm(mpg ~ am_f * vs_f + wt + hp, data = d)
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"),
              holds = list(wt = 3, hp = 150), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── glm ────────────────────────────────────────────────────────────────────
test_that("glm: factor*factor, x1 + x2", {
  fit <- glm(vs ~ am_f * cyl_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "cyl_f", x2vals = c("4", "6", "8"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("glm: factor*factor, x1 only", {
  fit <- glm(vs ~ am_f * cyl_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr ───────────────────────────────────────────────────────────────────
test_that("polr: factor*factor, x1 + x2", {
  fit <- polr(gear_o ~ am_f * vs_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor*factor, cut, x1 + x2", {
  fit <- polr(gear_o ~ am_f * vs_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"),
              cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr cloglog: factor*factor, cut, x1 + x2", {
  fit <- polr(gear_o ~ am_f * vs_f + wt, data = d, method = "cloglog")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "vs_f", x2vals = c("straight", "vee"),
              cut = 1, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: 3-level * 2-level factor, x1 + x2", {
  fit <- polr(gear_o ~ cyl_f * am_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "cyl_f", x1vals = c("4", "6", "8"),
              x2name = "am_f", x2vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("polr: factor*factor, x1 only (interaction in model)", {
  fit <- polr(gear_o ~ am_f * vs_f + wt, data = d, method = "logistic")
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})


# ============================================================================
# 6. SVYGLM
# ============================================================================

test_that("svyglm: no x (marginal)", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ wt + hp, design = des, family = binomial())
  res <- post(fit, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("svyglm: x1 only", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ wt + hp, design = des, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("svyglm: x1 + x2", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ wt + hp, design = des, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("svyglm: factor x1", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ am_f + wt, design = des, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("svyglm: with post-estimation weights argument", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ wt + hp, design = des, family = binomial())
  w <- rep(1, nrow(d))
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), weights = w, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("svyglm: factor*numeric interaction, x1 + x2", {
  skip_if_not_installed("survey")
  library(survey)
  d$svywt <- runif(nrow(d), 0.5, 2)
  des <- svydesign(ids = ~1, weights = ~svywt, data = d)
  fit <- svyglm(vs ~ am_f * wt, design = des, family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
              x2name = "wt", x2vals = c(2, 4), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})


# ============================================================================
# 7. BUG FIX VERIFICATION
# ============================================================================

# ── Bug 1: cloglog cut correctness ────────────────────────────────────────
test_that("Bug1: polr cloglog cut probabilities are in [0, 1]", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  set.seed(42)
  res <- post(fit, cut = 1, n.sims = N)
  expect_true(res@est[1, 1] >= 0 && res@est[1, 1] <= 1)
})

test_that("Bug1: polr cloglog cut x1 probabilities are in [0, 1]", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  set.seed(42)
  res <- post(fit, x1name = "wt", x1vals = c(2, 3, 4), cut = 1, n.sims = N)
  for (i in 1:3) {
    expect_true(res@est[i, 1] >= 0 && res@est[i, 1] <= 1,
                info = paste("row", i, "value:", res@est[i, 1]))
  }
})

test_that("Bug1: cloglog cut matches manual 1-link(tau-eta) computation", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  set.seed(99)
  res_cut   <- post(fit, cut = 1, n.sims = 1000, digits = 6)
  set.seed(99)
  res_nocut <- post(fit, n.sims = 1000, digits = 6)
  # P(Y > 1) should equal sum of P(Y=2) + P(Y=3) from the no-cut branch
  p_gt1 <- unname(res_cut@est[1, 1])
  p_y2_plus_y3 <- unname(res_nocut@est[2, 1] + res_nocut@est[3, 1])
  expect_equal(p_gt1, p_y2_plus_y3, tolerance = 0.01)
})

test_that("Bug1: logistic cut matches no-cut sum (sanity check)", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  set.seed(99)
  res_cut   <- post(fit, cut = 1, n.sims = 1000, digits = 6)
  set.seed(99)
  res_nocut <- post(fit, n.sims = 1000, digits = 6)
  p_gt1 <- unname(res_cut@est[1, 1])
  p_y2_plus_y3 <- unname(res_nocut@est[2, 1] + res_nocut@est[3, 1])
  expect_equal(p_gt1, p_y2_plus_y3, tolerance = 0.01)
})

# ── Bug 2: polr labels ───────────────────────────────────────────────────
test_that("Bug2: polr labels show actual levels, not 1:n", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, n.sims = N)
  labs <- rownames(res@est)
  # gear_o levels are "3","4","5"
  expect_true(all(grepl("3|4|5", labs)))
  expect_false(any(grepl("Y = 1$|Y = 2$", labs)))
})

test_that("Bug2: polr x1 branch labels show actual levels in 3rd dim", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = N)
  ylabs <- dimnames(res@est)[[3]]
  expect_true(all(grepl("3|4|5", ylabs)))
})

test_that("Bug2: polr x1+x2 branch labels show actual levels in 4th dim", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  ylabs <- dimnames(res@est)[[4]]
  expect_true(all(grepl("3|4|5", ylabs)))
})

test_that("Bug2: polr x1+x2 did labels show actual levels", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), n.sims = N)
  dlabs <- dimnames(res@did)[[1]]
  expect_true(all(grepl("3|4|5", dlabs)))
})

test_that("Bug2: polr cut x1+x2 did labels show actual levels", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  # cut branch x1+x2 did is a single row, labels should still be ok
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200), cut = 1, n.sims = N)
  expect_false(is.null(res@did))
  expect_false(any(is.na(res@did)))
})

# ── Bug 3: holds name validation ─────────────────────────────────────────
test_that("Bug3: misspelled holds name errors in glm", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "wt", x1vals = c(2, 4),
         holds = list(hpp = 100), n.sims = 20),
    "holds variable.*not found"
  )
})

test_that("Bug3: misspelled holds name errors in polr", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  expect_error(
    post(fit, holds = list(wtt = 3), n.sims = 20),
    "holds variable.*not found"
  )
})

test_that("Bug3: misspelled holds name errors in lm", {
  fit <- lm(mpg ~ wt + hp, data = d)
  expect_error(
    post(fit, holds = list(Weight = 3), n.sims = 20),
    "holds variable.*not found"
  )
})

test_that("Bug3: multiple holds with one bad name errors", {
  fit <- glm(vs ~ wt + hp + am, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "wt", x1vals = c(2, 4),
         holds = list(hp = 100, ammm = 1), n.sims = 20),
    "holds variable.*not found"
  )
})

# ── Bug 4: holds factor level validation ─────────────────────────────────
test_that("Bug4: invalid factor level in glm holds errors", {
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "wt", x1vals = c(2, 4),
         holds = list(am_f = "INVALID"), n.sims = 20),
    "holds value.*not a valid level"
  )
})

test_that("Bug4: invalid factor level in polr holds errors", {
  fit <- polr(gear_o ~ am_f + wt, data = d, method = "logistic")
  expect_error(
    post(fit, holds = list(am_f = "WRONG"), n.sims = 20),
    "holds value.*not a valid level"
  )
})

test_that("Bug4: valid factor level in holds works", {
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              holds = list(am_f = "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── Bug 5: single-quantile aperm ─────────────────────────────────────────
test_that("Bug5: single quantile in glm x1+x2 branch", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              quantiles = 0.5, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_equal(ncol(res@est), 2)  # mean + 1 quantile
})

test_that("Bug5: single quantile in glm x1 only", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              quantiles = 0.5, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

test_that("Bug5: single quantile in lm x1+x2 branch", {
  fit <- lm(mpg ~ wt + hp, data = d)
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              quantiles = 0.5, n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── Bug 6: sims slot array shape ─────────────────────────────────────────
test_that("Bug6: polr cut no-variable sims slot is 2D array", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, cut = 1, n.sims = N)
  expect_true(is.array(res@sims))
  expect_equal(length(dim(res@sims)), 2)
})


# ============================================================================
# 8. ADDITIONAL EDGE CASES AND STRESS TESTS
# ============================================================================

# ── weight validation ────────────────────────────────────────────────────
test_that("wrong-length glm weights error", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  expect_error(
    post(fit, weights = rep(1, 2), n.sims = 20),
    "weights must have the same length"
  )
})

test_that("wrong-length polr weights error", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  expect_error(
    post(fit, weights = rep(1, 3), n.sims = 20),
    "weights must have the same length"
  )
})

# ── variable name validation ─────────────────────────────────────────────
test_that("misspelled x1name in glm errors", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "weight", x1vals = c(2, 4), n.sims = 20),
    "x1name='weight' is not a variable"
  )
})

test_that("misspelled x2name in glm errors", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "wt", x1vals = c(2, 4),
         x2name = "horsepower", x2vals = c(100, 200), n.sims = 20),
    "x2name='horsepower' is not a variable"
  )
})

test_that("misspelled x1name in polr errors", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  expect_error(
    post(fit, x1name = "weight", x1vals = c(2, 4), n.sims = 20),
    "x1name='weight' is not a variable"
  )
})

# ── inline transformation detection ──────────────────────────────────────
test_that("inline factor(x) in x1name errors", {
  fit <- glm(vs ~ factor(am) + wt, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "factor(am)", x1vals = c("0", "1"), n.sims = 20),
    "appears to use an inline transformation"
  )
})

# ── invalid factor levels in x1vals/x2vals ───────────────────────────────
test_that("invalid factor x1vals errors with level list", {
  fit <- glm(vs ~ am_f + wt, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "am_f", x1vals = c("auto", "INVALID"), n.sims = 20),
    "x1vals must be valid levels"
  )
})

test_that("invalid factor x2vals errors with level list", {
  fit <- glm(vs ~ am_f + vs_f, data = d, family = binomial())
  expect_error(
    post(fit, x1name = "am_f", x1vals = c("auto", "manual"),
         x2name = "vs_f", x2vals = c("straight", "NOPE"), n.sims = 20),
    "x2vals must be valid levels"
  )
})

# ── numeric-like factor levels passed as numbers ─────────────────────────
test_that("numeric factor levels coerced correctly", {
  fit <- glm(vs ~ cyl_f + wt, data = d, family = binomial())
  res <- post(fit, x1name = "cyl_f", x1vals = c(4, 6, 8), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── three quantiles ──────────────────────────────────────────────────────
test_that("three quantiles work", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              quantiles = c(0.025, 0.5, 0.975), n.sims = N)
  expect_s4_class(res, "post")
  expect_equal(ncol(res@est), 4)  # mean + 3 quantiles
})

# ── many x1vals ──────────────────────────────────────────────────────────
test_that("many x1vals works", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  vals <- seq(2, 5, by = 0.5)
  res <- post(fit, x1name = "wt", x1vals = vals, n.sims = N)
  expect_equal(nrow(res@est), length(vals) + 1)
})

# ── three x2vals ─────────────────────────────────────────────────────────
test_that("three x2vals works in glm", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 150, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_equal(dim(res@est)[3], 3)
})

test_that("three x2vals works in polr", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 150, 200), n.sims = N)
  expect_s4_class(res, "post")
  expect_equal(dim(res@est)[3], 3)
})

# ── custom did parameter ─────────────────────────────────────────────────
test_that("custom did selects correct x2 pair", {
  fit <- glm(vs ~ wt + hp, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 150, 200),
              did = c(100, 200), n.sims = N)
  expect_false(is.null(res@did))
  expect_false(any(is.na(res@did)))
})

# ── multiple holds ───────────────────────────────────────────────────────
test_that("multiple holds work together", {
  fit <- glm(vs ~ wt + hp + am + drat, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              holds = list(am = 1, drat = 3.5), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── factor defined via transform() in model call ─────────────────────────
test_that("factor created via transform works", {
  fit <- glm(vs ~ am_f + wt,
             data = transform(mtcars, am_f = factor(ifelse(am == 1, "manual", "auto"))),
             family = binomial())
  res <- post(fit, x1name = "am_f", x1vals = c("auto", "manual"), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr with different ordered response ─────────────────────────────────
test_that("polr with cyl_o (levels 4,6,8) has correct labels", {
  fit <- polr(cyl_o ~ drat + qsec, data = d, method = "logistic")
  res <- post(fit, n.sims = N)
  labs <- rownames(res@est)
  expect_true(all(grepl("4|6|8", labs)))
  expect_false(any(grepl("Y = 1$|Y = 2$|Y = 3$", labs)))
})

# ── k not clobbered by holds in x1+x2 branches ──────────────────────────
test_that("holds loop does not shadow k in glm x1+x2", {
  fit <- glm(vs ~ wt + hp + am, data = d, family = binomial())
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(am = 1), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
  expect_false(is.null(res@did))
})

test_that("holds loop does not shadow k in polr x1+x2", {
  fit <- polr(gear_o ~ wt + hp + am, data = d, method = "logistic")
  res <- post(fit, x1name = "wt", x1vals = c(2, 4),
              x2name = "hp", x2vals = c(100, 200),
              holds = list(am = 1), n.sims = N)
  expect_s4_class(res, "post")
  expect_false(any(is.na(res@est)))
})

# ── polr probabilities sum to 1 ─────────────────────────────────────────
test_that("polr category probabilities sum to ~1", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "logistic")
  set.seed(42)
  res <- post(fit, n.sims = 500)
  prob_sum <- sum(res@est[, 1])
  expect_equal(prob_sum, 1.0, tolerance = 0.02)
})

test_that("polr cloglog category probabilities sum to ~1", {
  fit <- polr(gear_o ~ wt + hp, data = d, method = "cloglog")
  set.seed(42)
  res <- post(fit, n.sims = 500)
  prob_sum <- sum(res@est[, 1])
  expect_equal(prob_sum, 1.0, tolerance = 0.02)
})


cat("\nAll tests passed!\n")
