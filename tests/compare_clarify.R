#!/usr/bin/env Rscript
# ================================================================
# POST() vs CLARIFY COMPARISON TEST SUITE
# ================================================================
# Run from project root:  Rscript tests/compare_clarify.R
# Output:                 tests/compare_clarify_report.txt

cat("================================================================\n")
cat("POST() vs CLARIFY COMPARISON TEST SUITE\n")
cat("================================================================\n\n")

# ── Section 0: Setup & Configuration ─────────────────────────────

for (pkg in c("MASS", "clarify")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Required package '%s' is not installed.", pkg))
}

library(MASS)
library(clarify)

source("postSim.R")
source("post.R")

N_SIMS <- 5000
SEED   <- 42

# Divergence thresholds
THRESH_MEAN <- 0.15
THRESH_SE   <- 0.15
THRESH_CI   <- 0.35
THRESH_ABS  <- 0.005
THRESH_LM_DEFAULT <- 0.50

# Data setup
d <- mtcars
d$am_f   <- factor(ifelse(d$am == 1, "manual", "auto"))
d$vs_f   <- factor(ifelse(d$vs == 1, "straight", "vee"))
d$cyl_f  <- factor(d$cyl)
d$gear_o <- ordered(d$gear)
d$cyl_o  <- ordered(d$cyl)

# Accumulators
all_results  <- list()
test_count   <- 0L
pass_count   <- 0L
fail_count   <- 0L
error_count  <- 0L
bench_results <- list()

# Link-function map
link_fns <- list(
  identity = identity,
  logit    = plogis,
  probit   = pnorm,
  cloglog  = function(x) 1 - exp(-exp(x)),
  log      = exp,
  logistic = plogis
)

# ── Section 1: Helper Functions ──────────────────────────────────

# Compare two simulation vectors; return 4-row data.frame
compare_sims <- function(post_vec, clarify_vec, label, qty, wider = FALSE) {
  p_mean <- mean(post_vec);  c_mean <- mean(clarify_vec)
  p_sd   <- sd(post_vec);    c_sd   <- sd(clarify_vec)
  p_ci   <- quantile(post_vec,    c(0.025, 0.975))
  c_ci   <- quantile(clarify_vec, c(0.025, 0.975))

  pooled_se <- (p_sd + c_sd) / 2
  thresh_m  <- if (wider) THRESH_LM_DEFAULT else THRESH_MEAN

  d_mean <- abs(p_mean - c_mean)
  r_mean <- if (pooled_se > 1e-8) d_mean / pooled_se else 0
  pass_mean <- r_mean < thresh_m || d_mean < THRESH_ABS

  mean_se <- pooled_se
  d_se    <- abs(p_sd - c_sd)
  r_se    <- if (mean_se > 1e-8) d_se / mean_se else 0
  pass_se <- r_se < THRESH_SE || d_se < THRESH_ABS

  d_lo <- abs(p_ci[1] - c_ci[1])
  d_hi <- abs(p_ci[2] - c_ci[2])
  r_lo <- if (pooled_se > 1e-8) d_lo / pooled_se else 0
  r_hi <- if (pooled_se > 1e-8) d_hi / pooled_se else 0
  pass_lo <- r_lo < THRESH_CI || d_lo < THRESH_ABS
  pass_hi <- r_hi < THRESH_CI || d_hi < THRESH_ABS

  data.frame(
    label       = rep(label, 4),
    quantity    = rep(qty, 4),
    metric      = c("mean", "SE", "CI_lo", "CI_hi"),
    post_val    = round(c(p_mean, p_sd, unname(p_ci[1]), unname(p_ci[2])), 6),
    clarify_val = round(c(c_mean, c_sd, unname(c_ci[1]), unname(c_ci[2])), 6),
    diff_ratio  = round(c(r_mean, r_se, r_lo, r_hi), 4),
    pass        = c(pass_mean, pass_se, pass_lo, pass_hi),
    stringsAsFactors = FALSE
  )
}

record <- function(df) {
  all_results[[length(all_results) + 1L]] <<- df
  n <- nrow(df)
  test_count <<- test_count + n
  pass_count <<- pass_count + sum(df$pass)
  fail_count <<- fail_count + sum(!df$pass)
}

# ── Factory: GLM/LM observed-value average via coefs ─────────────
make_glm_ova_fn <- function(fit, x1name = NULL, x1vals = NULL,
                            x2name = NULL, x2vals = NULL,
                            holds = NULL, link_fn, weights = NULL) {
  mdata  <- data.frame(fit$model)
  mterms <- terms(fit)
  mxlev  <- fit$xlevels
  n_obs  <- nrow(mdata)
  wi     <- if (is.null(weights)) rep(1, n_obs) else weights
  if (!is.null(x1name) && is.factor(mdata[[x1name]])) x1vals <- as.character(x1vals)
  if (!is.null(x2name) && is.factor(mdata[[x2name]])) x2vals <- as.character(x2vals)

  function(coefs) {
    res <- c()
    pred_at <- function(nd) {
      X <- model.matrix(mterms, data = nd, xlev = mxlev)
      weighted.mean(link_fn(X %*% coefs), wi)
    }
    if (is.null(x1name)) {
      nd <- mdata
      if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
      res <- c(marginal = pred_at(nd))
    } else if (is.null(x2name)) {
      for (v in x1vals) {
        nd <- mdata
        if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
        nd[[x1name]] <- v
        res <- c(res, setNames(pred_at(nd), paste0(x1name, "=", v)))
      }
    } else {
      for (x2v in x2vals) {
        for (x1v in x1vals) {
          nd <- mdata
          if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
          nd[[x1name]] <- x1v
          nd[[x2name]] <- x2v
          nm <- paste0(x1name, "=", x1v, ",", x2name, "=", x2v)
          res <- c(res, setNames(pred_at(nd), nm))
        }
      }
    }
    res
  }
}

# ── Factory: polr all-category probabilities via coefs ───────────
make_polr_probs_fn <- function(fit, x1name = NULL, x1vals = NULL,
                               x2name = NULL, x2vals = NULL,
                               holds = NULL, link_fn, weights = NULL) {
  mdata      <- data.frame(fit$model)
  mterms     <- terms(fit)
  mxlev      <- fit$xlevels
  beta_names <- names(coef(fit))
  zeta_names <- names(fit$zeta)
  y_levels   <- levels(fit$model[[1]])
  n_y        <- length(y_levels)
  n_obs      <- nrow(mdata)
  wi         <- if (is.null(weights)) rep(1, n_obs) else weights
  if (!is.null(x1name) && is.factor(mdata[[x1name]])) x1vals <- as.character(x1vals)
  if (!is.null(x2name) && is.factor(mdata[[x2name]])) x2vals <- as.character(x2vals)

  function(coefs) {
    beta <- coefs[beta_names]
    zeta <- coefs[zeta_names]
    tau  <- c(-Inf, zeta, Inf)

    probs_at <- function(nd) {
      X <- model.matrix(mterms, data = nd, xlev = mxlev)[, -1, drop = FALSE]
      eta <- X %*% beta
      pr <- numeric(n_y)
      for (z in seq_len(n_y))
        pr[z] <- weighted.mean(link_fn(tau[z + 1] - eta) - link_fn(tau[z] - eta), wi)
      pr
    }

    res <- c()
    if (is.null(x1name)) {
      nd <- mdata
      if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
      pr <- probs_at(nd)
      for (z in seq_len(n_y))
        res <- c(res, setNames(pr[z], paste0("Y=", y_levels[z])))
    } else if (is.null(x2name)) {
      for (v in x1vals) {
        nd <- mdata
        if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
        nd[[x1name]] <- v
        pr <- probs_at(nd)
        for (z in seq_len(n_y))
          res <- c(res, setNames(pr[z], paste0(x1name, "=", v, ",Y=", y_levels[z])))
      }
    } else {
      for (x2v in x2vals) {
        for (x1v in x1vals) {
          nd <- mdata
          if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
          nd[[x1name]] <- x1v
          nd[[x2name]] <- x2v
          pr <- probs_at(nd)
          for (z in seq_len(n_y))
            res <- c(res, setNames(pr[z],
              paste0(x1name, "=", x1v, ",", x2name, "=", x2v, ",Y=", y_levels[z])))
        }
      }
    }
    res
  }
}

# ── Factory: polr P(Y > cut) via coefs ──────────────────────────
make_polr_cut_fn <- function(fit, x1name = NULL, x1vals = NULL,
                             x2name = NULL, x2vals = NULL,
                             holds = NULL, link_fn, cut,
                             weights = NULL) {
  mdata      <- data.frame(fit$model)
  mterms     <- terms(fit)
  mxlev      <- fit$xlevels
  beta_names <- names(coef(fit))
  zeta_names <- names(fit$zeta)
  n_obs      <- nrow(mdata)
  wi         <- if (is.null(weights)) rep(1, n_obs) else weights
  if (!is.null(x1name) && is.factor(mdata[[x1name]])) x1vals <- as.character(x1vals)
  if (!is.null(x2name) && is.factor(mdata[[x2name]])) x2vals <- as.character(x2vals)

  function(coefs) {
    beta <- coefs[beta_names]
    zeta <- coefs[zeta_names]
    tau  <- c(-Inf, zeta, Inf)

    pgt_at <- function(nd) {
      X <- model.matrix(mterms, data = nd, xlev = mxlev)[, -1, drop = FALSE]
      eta <- X %*% beta
      weighted.mean(1 - link_fn(tau[cut + 1] - eta), wi)
    }

    res <- c()
    if (is.null(x1name)) {
      nd <- mdata
      if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
      res <- c("P(Y>cut)" = pgt_at(nd))
    } else if (is.null(x2name)) {
      for (v in x1vals) {
        nd <- mdata
        if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
        nd[[x1name]] <- v
        res <- c(res, setNames(pgt_at(nd), paste0(x1name, "=", v)))
      }
    } else {
      for (x2v in x2vals) {
        for (x1v in x1vals) {
          nd <- mdata
          if (!is.null(holds)) for (h in names(holds)) nd[[h]] <- holds[[h]]
          nd[[x1name]] <- x1v
          nd[[x2name]] <- x2v
          nm <- paste0(x1name, "=", x1v, ",", x2name, "=", x2v)
          res <- c(res, setNames(pgt_at(nd), nm))
        }
      }
    }
    res
  }
}

# ── Unified comparison runner ────────────────────────────────────
run_comparison <- function(id, desc, fit,
                           type     = c("glm", "polr"),
                           x1name   = NULL, x1vals = NULL,
                           x2name   = NULL, x2vals = NULL,
                           holds    = NULL, cut = NULL,
                           lm_normal = TRUE, wider = FALSE,
                           weights  = NULL) {
  type <- match.arg(type)
  cat(sprintf("  %-4s %-45s ", id, desc))

  tryCatch({
    # ── post() ──
    set.seed(SEED)
    p <- post(fit, x1name = x1name, x1vals = x1vals,
              x2name = x2name, x2vals = x2vals,
              holds = holds, cut = cut, weights = weights,
              n.sims = N_SIMS, digits = 10)

    # ── clarify sim() ──
    set.seed(SEED)
    if (type == "polr") {
      polr_coefs <- c(fit$coefficients, fit$zeta)
      s <- sim(fit, n = N_SIMS, coefs = polr_coefs, vcov = vcov(fit))
    } else {
      is_lm <- inherits(fit, "lm") && !inherits(fit, "glm")
      if (is_lm && lm_normal) {
        s <- sim(fit, n = N_SIMS, dist = "normal")
      } else {
        s <- sim(fit, n = N_SIMS)
      }
    }

    # ── clarify sim_apply ──
    if (type == "glm") {
      lname <- family(fit)$link
      lf    <- link_fns[[lname]]
      fn    <- make_glm_ova_fn(fit, x1name, x1vals, x2name, x2vals,
                               holds, lf, weights)
    } else if (is.null(cut)) {
      lf <- link_fns[[fit$method]]
      fn <- make_polr_probs_fn(fit, x1name, x1vals, x2name, x2vals,
                               holds, lf, weights)
    } else {
      lf <- link_fns[[fit$method]]
      fn <- make_polr_cut_fn(fit, x1name, x1vals, x2name, x2vals,
                             holds, lf, cut, weights)
    }
    c_res <- sim_apply(s, fn, verbose = FALSE)
    c_mat <- as.matrix(c_res)

    # ── Extract & compare ──
    results <- list()
    y_levels <- if (type == "polr") levels(fit$model[[1]]) else NULL

    if (type == "glm") {
      if (is.null(x1name)) {
        results[[1]] <- compare_sims(p@sims[, 1], c_mat[, 1],
                                     id, "marginal", wider)
      } else if (is.null(x2name)) {
        fv <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        for (i in seq_along(x1vals)) {
          cname <- paste0(x1name, "=", fv[i])
          results[[length(results) + 1L]] <- compare_sims(
            p@sims[, i], c_mat[, cname], id, cname, wider)
        }
      } else {
        fv1 <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        fv2 <- if (is.factor(fit$model[[x2name]])) as.character(x2vals) else x2vals
        for (j in seq_along(x2vals)) {
          for (i in seq_along(x1vals)) {
            cname <- paste0(x1name, "=", fv1[i], ",", x2name, "=", fv2[j])
            results[[length(results) + 1L]] <- compare_sims(
              p@sims[, i, j], c_mat[, cname], id, cname, wider)
          }
        }
      }

    } else if (is.null(cut)) {
      # polr all probs
      n_y <- length(y_levels)
      if (is.null(x1name)) {
        for (z in seq_len(n_y)) {
          cname <- paste0("Y=", y_levels[z])
          results[[length(results) + 1L]] <- compare_sims(
            p@sims[, z], c_mat[, cname], id, cname, wider)
        }
      } else if (is.null(x2name)) {
        fv <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        for (i in seq_along(x1vals)) {
          for (z in seq_len(n_y)) {
            cname <- paste0(x1name, "=", fv[i], ",Y=", y_levels[z])
            results[[length(results) + 1L]] <- compare_sims(
              p@sims[, i, z], c_mat[, cname], id, cname, wider)
          }
        }
      } else {
        fv1 <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        fv2 <- if (is.factor(fit$model[[x2name]])) as.character(x2vals) else x2vals
        for (j in seq_along(x2vals)) {
          for (i in seq_along(x1vals)) {
            for (z in seq_len(n_y)) {
              cname <- paste0(x1name, "=", fv1[i], ",", x2name, "=", fv2[j],
                              ",Y=", y_levels[z])
              results[[length(results) + 1L]] <- compare_sims(
                p@sims[, i, j, z], c_mat[, cname], id, cname, wider)
            }
          }
        }
      }

    } else {
      # polr P(Y > cut)
      if (is.null(x1name)) {
        results[[1]] <- compare_sims(p@sims[, 1], c_mat[, 1],
                                     id, "P(Y>cut)", wider)
      } else if (is.null(x2name)) {
        fv <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        for (i in seq_along(x1vals)) {
          cname <- paste0(x1name, "=", fv[i])
          results[[length(results) + 1L]] <- compare_sims(
            p@sims[, i], c_mat[, cname], id, cname, wider)
        }
      } else {
        fv1 <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals
        fv2 <- if (is.factor(fit$model[[x2name]])) as.character(x2vals) else x2vals
        for (j in seq_along(x2vals)) {
          for (i in seq_along(x1vals)) {
            cname <- paste0(x1name, "=", fv1[i], ",", x2name, "=", fv2[j])
            results[[length(results) + 1L]] <- compare_sims(
              p@sims[, i, j], c_mat[, cname], id, cname, wider)
          }
        }
      }
    }

    df <- do.call(rbind, results)
    record(df)
    n_fail <- sum(!df$pass)
    if (n_fail == 0) cat("PASS\n") else cat(sprintf("FAIL (%d/%d)\n", n_fail, nrow(df)))

  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
  })
}

# ── First-difference / DiD comparison runner ─────────────────────
run_fd_comparison <- function(id, desc, fit,
                              type     = c("glm", "polr"),
                              x1name, x1vals,
                              x2name   = NULL, x2vals = NULL,
                              holds    = NULL, cut = NULL,
                              lm_normal = TRUE, wider = FALSE) {
  type <- match.arg(type)
  cat(sprintf("  %-4s %-45s ", id, desc))

  tryCatch({
    set.seed(SEED)
    p <- post(fit, x1name = x1name, x1vals = x1vals,
              x2name = x2name, x2vals = x2vals,
              holds = holds, cut = cut,
              n.sims = N_SIMS, digits = 10)

    set.seed(SEED)
    if (type == "polr") {
      s <- sim(fit, n = N_SIMS,
               coefs = c(fit$coefficients, fit$zeta),
               vcov  = vcov(fit))
    } else {
      is_lm <- inherits(fit, "lm") && !inherits(fit, "glm")
      if (is_lm && lm_normal) s <- sim(fit, n = N_SIMS, dist = "normal")
      else                    s <- sim(fit, n = N_SIMS)
    }

    # Build sim_apply function (same as run_comparison)
    if (type == "glm") {
      lf <- link_fns[[family(fit)$link]]
      fn <- make_glm_ova_fn(fit, x1name, x1vals, x2name, x2vals, holds, lf)
    } else if (is.null(cut)) {
      lf <- link_fns[[fit$method]]
      fn <- make_polr_probs_fn(fit, x1name, x1vals, x2name, x2vals, holds, lf)
    } else {
      lf <- link_fns[[fit$method]]
      fn <- make_polr_cut_fn(fit, x1name, x1vals, x2name, x2vals, holds, lf, cut)
    }
    c_res <- sim_apply(s, fn, verbose = FALSE)
    c_mat <- as.matrix(c_res)

    results <- list()
    n_x1 <- length(x1vals)
    fv1 <- if (is.factor(fit$model[[x1name]])) as.character(x1vals) else x1vals

    if (is.null(x2name)) {
      # First difference only
      if (type == "glm" || !is.null(cut)) {
        # post FD sims
        post_fd <- p@sims[, n_x1] - p@sims[, 1]
        # clarify FD sims
        cn_hi <- paste0(x1name, "=", fv1[n_x1])
        cn_lo <- paste0(x1name, "=", fv1[1])
        clar_fd <- c_mat[, cn_hi] - c_mat[, cn_lo]
        results[[1]] <- compare_sims(post_fd, clar_fd, id, "FD", wider)
      } else {
        # polr all probs: FD per category
        y_levels <- levels(fit$model[[1]])
        for (z in seq_along(y_levels)) {
          post_fd <- p@sims[, n_x1, z] - p@sims[, 1, z]
          cn_hi <- paste0(x1name, "=", fv1[n_x1], ",Y=", y_levels[z])
          cn_lo <- paste0(x1name, "=", fv1[1],    ",Y=", y_levels[z])
          clar_fd <- c_mat[, cn_hi] - c_mat[, cn_lo]
          results[[length(results) + 1L]] <- compare_sims(
            post_fd, clar_fd, id, paste0("FD,Y=", y_levels[z]), wider)
        }
      }

    } else {
      # DiD
      fv2 <- if (is.factor(fit$model[[x2name]])) as.character(x2vals) else x2vals
      n_x2 <- length(x2vals)

      if (type == "glm" || !is.null(cut)) {
        # post DiD sims
        post_did <- (p@sims[, n_x1, n_x2] - p@sims[, 1, n_x2]) -
                    (p@sims[, n_x1, 1]    - p@sims[, 1, 1])
        # clarify DiD sims
        cn_hh <- paste0(x1name, "=", fv1[n_x1], ",", x2name, "=", fv2[n_x2])
        cn_lh <- paste0(x1name, "=", fv1[1],    ",", x2name, "=", fv2[n_x2])
        cn_hl <- paste0(x1name, "=", fv1[n_x1], ",", x2name, "=", fv2[1])
        cn_ll <- paste0(x1name, "=", fv1[1],    ",", x2name, "=", fv2[1])
        clar_did <- (c_mat[, cn_hh] - c_mat[, cn_lh]) -
                    (c_mat[, cn_hl] - c_mat[, cn_ll])
        results[[1]] <- compare_sims(post_did, clar_did, id, "DiD", wider)
      } else {
        # polr all probs: DiD per category
        y_levels <- levels(fit$model[[1]])
        for (z in seq_along(y_levels)) {
          post_did <- (p@sims[, n_x1, n_x2, z] - p@sims[, 1, n_x2, z]) -
                      (p@sims[, n_x1, 1, z]    - p@sims[, 1, 1, z])
          cn_hh <- paste0(x1name,"=",fv1[n_x1],",",x2name,"=",fv2[n_x2],",Y=",y_levels[z])
          cn_lh <- paste0(x1name,"=",fv1[1],   ",",x2name,"=",fv2[n_x2],",Y=",y_levels[z])
          cn_hl <- paste0(x1name,"=",fv1[n_x1],",",x2name,"=",fv2[1],   ",Y=",y_levels[z])
          cn_ll <- paste0(x1name,"=",fv1[1],   ",",x2name,"=",fv2[1],   ",Y=",y_levels[z])
          clar_did <- (c_mat[, cn_hh] - c_mat[, cn_lh]) -
                      (c_mat[, cn_hl] - c_mat[, cn_ll])
          results[[length(results) + 1L]] <- compare_sims(
            post_did, clar_did, id, paste0("DiD,Y=", y_levels[z]), wider)
        }
      }
    }

    df <- do.call(rbind, results)
    record(df)
    n_fail <- sum(!df$pass)
    if (n_fail == 0) cat("PASS\n") else cat(sprintf("FAIL (%d/%d)\n", n_fail, nrow(df)))

  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
  })
}


# ================================================================
# Section 2: Test Groups A – H
# ================================================================

# ── Group A: No factors, no interactions ─────────────────────────
cat("\nGroup A: No factors, no interactions\n")

run_comparison("A1", "lm, no x, marginal",
  fit = lm(mpg ~ wt + hp, data = d), type = "glm")

run_comparison("A2", "lm, x1=wt",
  fit = lm(mpg ~ wt + hp, data = d), type = "glm",
  x1name = "wt", x1vals = c(2, 4))

run_comparison("A3", "lm, x1+x2 (wt, hp)",
  fit = lm(mpg ~ wt + hp, data = d), type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))

run_comparison("A4", "glm logit, no x",
  fit = glm(vs ~ wt + hp, data = d, family = binomial()), type = "glm")

run_comparison("A5", "glm logit, x1=wt",
  fit = glm(vs ~ wt + hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 3, 4))

run_comparison("A6", "glm logit, x1+x2 (wt, hp)",
  fit = glm(vs ~ wt + hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))

run_comparison("A7", "glm probit, x1=wt",
  fit = glm(vs ~ wt + hp, data = d, family = binomial(link = "probit")),
  type = "glm", x1name = "wt", x1vals = c(2, 4))

run_comparison("A8", "glm cloglog, x1=wt",
  fit = glm(vs ~ wt + hp, data = d, family = binomial(link = "cloglog")),
  type = "glm", x1name = "wt", x1vals = c(2, 4))

run_comparison("A9", "glm poisson, x1=wt",
  fit = glm(carb ~ wt + hp, data = d, family = poisson()),
  type = "glm", x1name = "wt", x1vals = c(2, 4))

run_comparison("A10", "polr logistic, no x, all probs",
  fit = polr(gear_o ~ wt + hp, data = d, method = "logistic"),
  type = "polr")


# ── Group B: Factor variables, no interactions ───────────────────
cat("\nGroup B: Factor variables, no interactions\n")

run_comparison("B1", "lm, factor x1=am_f",
  fit = lm(mpg ~ am_f + wt, data = d), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"))

run_comparison("B2", "glm logit, factor x1=am_f",
  fit = glm(vs ~ am_f + wt, data = d, family = binomial()), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"))

run_comparison("B3", "glm logit, factor x1=cyl_f (3 levels)",
  fit = glm(vs ~ cyl_f + wt, data = d, family = binomial()), type = "glm",
  x1name = "cyl_f", x1vals = c("4", "6", "8"))

run_comparison("B4", "lm, factor x1 + numeric x2",
  fit = lm(mpg ~ am_f + wt, data = d), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "wt", x2vals = c(2, 4))

run_comparison("B5", "polr logistic, factor x1, all probs",
  fit = polr(gear_o ~ am_f + wt, data = d, method = "logistic"),
  type = "polr",
  x1name = "am_f", x1vals = c("auto", "manual"))

run_comparison("B6", "polr logistic, factor x1, cut=1",
  fit = polr(gear_o ~ am_f + wt, data = d, method = "logistic"),
  type = "polr", cut = 1,
  x1name = "am_f", x1vals = c("auto", "manual"))


# ── Group C: Numeric interactions ────────────────────────────────
cat("\nGroup C: Numeric interactions\n")

run_comparison("C1", "lm, wt*hp, x1=wt",
  fit = lm(mpg ~ wt * hp, data = d), type = "glm",
  x1name = "wt", x1vals = c(2, 4))

run_comparison("C2", "glm logit, wt*hp, x1=wt",
  fit = glm(vs ~ wt * hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4))

run_comparison("C3", "glm logit, wt*hp, x1+x2",
  fit = glm(vs ~ wt * hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))

run_comparison("C4", "polr logistic, wt*am, x1=wt, all probs",
  fit = polr(gear_o ~ wt * am, data = d, method = "logistic"),
  type = "polr",
  x1name = "wt", x1vals = c(2, 4))


# ── Group D: Factor x numeric interactions ───────────────────────
cat("\nGroup D: Factor x numeric interactions\n")

run_comparison("D1", "lm, am_f*wt, factor x1 + numeric x2",
  fit = lm(mpg ~ am_f * wt, data = d), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "wt", x2vals = c(2, 4))

run_comparison("D2", "glm logit, am_f*wt, factor x1 + num x2",
  fit = glm(vs ~ am_f * wt, data = d, family = binomial()), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "wt", x2vals = c(2, 4))

run_comparison("D3", "polr logistic, am_f*wt, fac x1 + num x2",
  fit = polr(gear_o ~ am_f * wt, data = d, method = "logistic"),
  type = "polr",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "wt", x2vals = c(2, 4))

run_comparison("D4", "glm logit, am_f*wt, num x1 + fac x2",
  fit = glm(vs ~ am_f * wt, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  x2name = "am_f", x2vals = c("auto", "manual"))


# ── Group E: Factor x factor interactions ────────────────────────
cat("\nGroup E: Factor x factor interactions\n")

run_comparison("E1", "lm, am_f*vs_f, factor x1 + factor x2",
  fit = lm(mpg ~ am_f * vs_f + wt, data = d), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "vs_f", x2vals = c("straight", "vee"))

run_comparison("E2", "glm logit, am_f*cyl_f, fac x1 + fac x2",
  fit = glm(vs ~ am_f * cyl_f + wt, data = d, family = binomial()),
  type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "cyl_f", x2vals = c("4", "6", "8"))

run_comparison("E3", "polr logistic, am_f*vs_f, fac x1+x2, cut=1",
  fit = polr(gear_o ~ am_f * vs_f + wt, data = d, method = "logistic"),
  type = "polr", cut = 1,
  x1name = "am_f", x1vals = c("auto", "manual"),
  x2name = "vs_f", x2vals = c("straight", "vee"))


# ── Group F: polr link functions ─────────────────────────────────
cat("\nGroup F: polr link functions\n")

run_comparison("F1", "polr probit, x1=wt, all probs",
  fit = polr(gear_o ~ wt + hp, data = d, method = "probit"),
  type = "polr",
  x1name = "wt", x1vals = c(2, 4))

run_comparison("F2", "polr cloglog, x1=wt, all probs",
  fit = polr(gear_o ~ wt + hp, data = d, method = "cloglog"),
  type = "polr",
  x1name = "wt", x1vals = c(2, 4))

run_comparison("F3", "polr probit, x1+x2, cut=1",
  fit = polr(gear_o ~ wt + hp, data = d, method = "probit"),
  type = "polr", cut = 1,
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))


# ── Group G: First differences & DiD ────────────────────────────
cat("\nGroup G: First differences & DiD\n")

run_fd_comparison("G1", "FD, glm logit, x1 only",
  fit = glm(vs ~ wt + hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4))

run_fd_comparison("G2", "FD, lm, factor x1",
  fit = lm(mpg ~ am_f + wt, data = d), type = "glm",
  x1name = "am_f", x1vals = c("auto", "manual"))

run_fd_comparison("G3", "DiD, glm logit, x1+x2",
  fit = glm(vs ~ wt + hp, data = d, family = binomial()), type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))

run_fd_comparison("G4", "DiD, polr logistic, cut=1",
  fit = polr(gear_o ~ wt + hp, data = d, method = "logistic"),
  type = "polr", cut = 1,
  x1name = "wt", x1vals = c(2, 4),
  x2name = "hp", x2vals = c(100, 200))


# ── Group H: Special cases ──────────────────────────────────────
cat("\nGroup H: Special cases\n")

# H1: holds parameter
run_comparison("H1", "glm logit, with holds",
  fit = glm(vs ~ wt + hp + drat, data = d, family = binomial()),
  type = "glm",
  x1name = "wt", x1vals = c(2, 4),
  holds = list(drat = 3.5))

# H2: custom weights
{
  cat(sprintf("  %-4s %-45s ", "H2", "glm logit, custom weights"))
  tryCatch({
    fit_h2 <- glm(vs ~ wt + hp, data = d, family = binomial())
    w_h2   <- runif(nrow(d), 0.5, 2)

    set.seed(SEED)
    p_h2 <- post(fit_h2, weights = w_h2, n.sims = N_SIMS, digits = 10)

    set.seed(SEED)
    s_h2 <- sim(fit_h2, n = N_SIMS)
    fn_h2 <- make_glm_ova_fn(fit_h2, link_fn = plogis, weights = w_h2)
    c_h2  <- sim_apply(s_h2, fn_h2, verbose = FALSE)
    cm_h2 <- as.matrix(c_h2)

    df_h2 <- compare_sims(p_h2@sims[, 1], cm_h2[, 1], "H2", "weighted_marginal")
    record(df_h2)
    if (all(df_h2$pass)) cat("PASS\n") else cat(sprintf("FAIL (%d/%d)\n",
      sum(!df_h2$pass), nrow(df_h2)))
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
  })
}

# H3: polr probs sum to 1
{
  cat(sprintf("  %-4s %-45s ", "H3", "polr probs sum to 1 (both)"))
  tryCatch({
    fit_h3 <- polr(gear_o ~ wt + hp, data = d, method = "logistic")

    set.seed(SEED)
    p_h3 <- post(fit_h3, n.sims = N_SIMS, digits = 10)

    set.seed(SEED)
    s_h3 <- sim(fit_h3, n = N_SIMS,
                coefs = c(fit_h3$coefficients, fit_h3$zeta),
                vcov  = vcov(fit_h3))
    fn_h3 <- make_polr_probs_fn(fit_h3, link_fn = plogis)
    c_h3  <- sim_apply(s_h3, fn_h3, verbose = FALSE)
    cm_h3 <- as.matrix(c_h3)

    # Check that probabilities sum to ~1 for every simulation
    post_sums    <- rowSums(p_h3@sims)
    clarify_sums <- rowSums(cm_h3)
    post_ok    <- all(abs(post_sums - 1) < 0.01)
    clarify_ok <- all(abs(clarify_sums - 1) < 0.01)

    df_h3 <- data.frame(
      label       = rep("H3", 2),
      quantity    = c("post_prob_sum", "clarify_prob_sum"),
      metric      = c("sum_to_1", "sum_to_1"),
      post_val    = round(c(mean(post_sums), mean(clarify_sums)), 6),
      clarify_val = round(c(1, 1), 6),
      diff_ratio  = round(c(abs(mean(post_sums) - 1), abs(mean(clarify_sums) - 1)), 6),
      pass        = c(post_ok, clarify_ok),
      stringsAsFactors = FALSE
    )
    record(df_h3)
    if (all(df_h3$pass)) cat("PASS\n") else cat(sprintf("FAIL (%d/%d)\n",
      sum(!df_h3$pass), nrow(df_h3)))
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
  })
}

# H4: lm with default dist (expected divergence — informational only)
# postSim uses chi-squared sigma + conditional MVN; clarify default uses
# multivariate t.  These differ by design so H4 is scored separately.
h4_results <- NULL
{
  cat(sprintf("  %-4s %-45s ", "H4", "lm default dist (expected divergence)"))
  .h4_ok <- tryCatch({
    fit_h4 <- lm(mpg ~ wt + hp, data = d)

    set.seed(SEED)
    p_h4 <- post(fit_h4, x1name = "wt", x1vals = c(2, 4),
                 n.sims = N_SIMS, digits = 10)

    set.seed(SEED)
    # Clarify with DEFAULT dist (t-distribution) instead of "normal"
    s_h4 <- sim(fit_h4, n = N_SIMS)
    fn_h4 <- make_glm_ova_fn(fit_h4, x1name = "wt", x1vals = c(2, 4),
                              link_fn = identity)
    c_h4  <- sim_apply(s_h4, fn_h4, verbose = FALSE)
    cm_h4 <- as.matrix(c_h4)

    tmp <- list()
    for (i in 1:2) {
      cname <- paste0("wt=", c(2, 4)[i])
      tmp[[i]] <- compare_sims(
        p_h4@sims[, i], cm_h4[, cname], "H4*", cname, wider = TRUE)
    }
    tmp   # return value
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
    NULL
  })
  if (!is.null(.h4_ok)) {
    h4_results <- do.call(rbind, .h4_ok)
    cat("INFO (expected divergence)\n")
  }
}

# H5: lm with dist="t" in post() vs clarify default t-distribution
# postSim dist="t" uses the standard MVt: mu + Z/sqrt(W).  clarify's rmvt
# divides (mu + Z) / sqrt(W), scaling the mean by the chi-squared factor.
# This non-standard implementation inflates variance proportional to mu^2,
# so H5 is scored separately like H4.
h5_results <- NULL
{
  cat(sprintf("  %-4s %-45s ", "H5", "lm dist='t' vs clarify default t"))
  .h5_ok <- tryCatch({
    fit_h5 <- lm(mpg ~ wt + hp, data = d)

    set.seed(SEED)
    p_h5 <- post(fit_h5, x1name = "wt", x1vals = c(2, 4),
                 n.sims = N_SIMS, digits = 10, dist = "t")

    set.seed(SEED)
    # Clarify with DEFAULT dist (t-distribution)
    s_h5 <- sim(fit_h5, n = N_SIMS)
    fn_h5 <- make_glm_ova_fn(fit_h5, x1name = "wt", x1vals = c(2, 4),
                              link_fn = identity)
    c_h5  <- sim_apply(s_h5, fn_h5, verbose = FALSE)
    cm_h5 <- as.matrix(c_h5)

    tmp <- list()
    for (i in 1:2) {
      cname <- paste0("wt=", c(2, 4)[i])
      tmp[[i]] <- compare_sims(
        p_h5@sims[, i], cm_h5[, cname], "H5*", cname, wider = TRUE)
    }
    tmp
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    error_count <<- error_count + 1L
    NULL
  })
  if (!is.null(.h5_ok)) {
    h5_results <- do.call(rbind, .h5_ok)
    cat("INFO (expected divergence)\n")
  }
}


# ================================================================
# Section 3: Benchmark Timing
# ================================================================
cat("\nBenchmarks\n")

N_BENCH  <- 1000
N_REPS   <- 3

bench_specs <- list(
  list(label = "lm, x1 only",
       fit_expr = quote(lm(mpg ~ wt + hp, data = d)),
       type = "glm", x1name = "wt", x1vals = c(2, 4),
       x2name = NULL, x2vals = NULL, cut = NULL),
  list(label = "glm logit, x1+x2",
       fit_expr = quote(glm(vs ~ wt + hp, data = d, family = binomial())),
       type = "glm", x1name = "wt", x1vals = c(2, 4),
       x2name = "hp", x2vals = c(100, 200), cut = NULL),
  list(label = "polr logistic, x1, all probs",
       fit_expr = quote(polr(gear_o ~ wt + hp, data = d, method = "logistic")),
       type = "polr", x1name = "wt", x1vals = c(2, 4),
       x2name = NULL, x2vals = NULL, cut = NULL),
  list(label = "polr logistic, x1+x2, cut=1",
       fit_expr = quote(polr(gear_o ~ wt + hp, data = d, method = "logistic")),
       type = "polr", x1name = "wt", x1vals = c(2, 4),
       x2name = "hp", x2vals = c(100, 200), cut = 1),
  list(label = "glm logit, fac*num interaction",
       fit_expr = quote(glm(vs ~ am_f * wt, data = d, family = binomial())),
       type = "glm", x1name = "am_f", x1vals = c("auto", "manual"),
       x2name = "wt", x2vals = c(2, 4), cut = NULL)
)

for (bs in bench_specs) {
  fit_b <- eval(bs$fit_expr)
  cat(sprintf("  %-40s ", bs$label))

  post_times    <- numeric(N_REPS)
  clarify_times <- numeric(N_REPS)

  for (r in seq_len(N_REPS)) {
    # Time post()
    set.seed(SEED + r)
    t0 <- proc.time()["elapsed"]
    suppressMessages(
      post(fit_b, x1name = bs$x1name, x1vals = bs$x1vals,
           x2name = bs$x2name, x2vals = bs$x2vals,
           cut = bs$cut, n.sims = N_BENCH, digits = 10)
    )
    post_times[r] <- proc.time()["elapsed"] - t0

    # Time clarify pipeline
    set.seed(SEED + r)
    t0 <- proc.time()["elapsed"]
    if (bs$type == "polr") {
      s_b <- sim(fit_b, n = N_BENCH,
                 coefs = c(fit_b$coefficients, fit_b$zeta),
                 vcov  = vcov(fit_b))
      lf_b <- link_fns[[fit_b$method]]
      if (is.null(bs$cut)) {
        fn_b <- make_polr_probs_fn(fit_b, bs$x1name, bs$x1vals,
                                   bs$x2name, bs$x2vals, link_fn = lf_b)
      } else {
        fn_b <- make_polr_cut_fn(fit_b, bs$x1name, bs$x1vals,
                                 bs$x2name, bs$x2vals, link_fn = lf_b,
                                 cut = bs$cut)
      }
    } else {
      s_b  <- sim(fit_b, n = N_BENCH)
      lf_b <- link_fns[[family(fit_b)$link]]
      fn_b <- make_glm_ova_fn(fit_b, bs$x1name, bs$x1vals,
                               bs$x2name, bs$x2vals, link_fn = lf_b)
    }
    sim_apply(s_b, fn_b, verbose = FALSE)
    clarify_times[r] <- proc.time()["elapsed"] - t0
  }

  med_post    <- median(post_times)
  med_clarify <- median(clarify_times)
  ratio       <- med_post / med_clarify

  bench_results[[length(bench_results) + 1L]] <- data.frame(
    model   = bs$label,
    post_s  = round(med_post, 3),
    clarify_s = round(med_clarify, 3),
    ratio   = round(ratio, 2),
    stringsAsFactors = FALSE
  )

  cat(sprintf("post=%.3fs  clarify=%.3fs  ratio=%.2f\n",
              med_post, med_clarify, ratio))
}


# ================================================================
# Section 4: Generate Report
# ================================================================
cat("\nGenerating report...\n")

report_path <- "tests/compare_clarify_report.txt"

all_df <- if (length(all_results) > 0) do.call(rbind, all_results) else
  data.frame(label = character(0), quantity = character(0),
             metric = character(0), post_val = numeric(0),
             clarify_val = numeric(0), diff_ratio = numeric(0),
             pass = logical(0), stringsAsFactors = FALSE)

bench_df <- if (length(bench_results) > 0) do.call(rbind, bench_results) else
  data.frame(model = character(0), post_s = numeric(0),
             clarify_s = numeric(0), ratio = numeric(0),
             stringsAsFactors = FALSE)

lines <- character(0)
add <- function(...) lines <<- c(lines, sprintf(...))

add("================================================================")
add("POST() vs CLARIFY COMPARISON REPORT")
add("================================================================")
add("Date:          %s", Sys.time())
add("n.sims:        %d", N_SIMS)
add("R version:     %s", R.version.string)
add("clarify:       %s", as.character(packageVersion("clarify")))
add("MASS:          %s", as.character(packageVersion("MASS")))
add("")
add("Thresholds:    mean < %.2f | SE < %.2f | CI < %.2f | abs < %.4f",
    THRESH_MEAN, THRESH_SE, THRESH_CI, THRESH_ABS)
add("               lm default dist (wider): mean < %.2f", THRESH_LM_DEFAULT)
add("")

# Accuracy section
add("ACCURACY COMPARISONS")
add("--------------------------------------------------------------")
add("%-5s %-38s %-6s %10s %10s %8s %s",
    "ID", "Quantity", "Metric", "post", "clarify", "d/SE", "Result")
add("--------------------------------------------------------------")

for (i in seq_len(nrow(all_df))) {
  r <- all_df[i, ]
  tag <- if (r$pass) "PASS" else "FAIL"
  add("%-5s %-38s %-6s %10.6f %10.6f %8.4f %s",
      r$label, r$quantity, r$metric, r$post_val, r$clarify_val,
      r$diff_ratio, tag)
}

add("")

# H4/H5 expected divergence section
h4_df <- h4_results
h5_df <- h5_results
div_df <- rbind(h4_df, h5_df)
if (!is.null(div_df) && nrow(div_df) > 0) {
  add("EXPECTED DIVERGENCE: H4 & H5 (lm t-distribution differences)")
  add("H4: postSim default (chi-sq/MVN) vs clarify default (non-standard MVt)")
  add("H5: postSim dist='t' (standard MVt) vs clarify default (non-standard MVt)")
  add("clarify's rmvt divides (mu+Z)/sqrt(W) instead of mu+Z/sqrt(W),")
  add("inflating variance proportional to mu^2.")
  add("--------------------------------------------------------------")
  add("%-5s %-38s %-6s %10s %10s %8s",
      "ID", "Quantity", "Metric", "post", "clarify", "d/SE")
  add("--------------------------------------------------------------")
  for (i in seq_len(nrow(div_df))) {
    r <- div_df[i, ]
    add("%-5s %-38s %-6s %10.6f %10.6f %8.4f",
        r$label, r$quantity, r$metric, r$post_val, r$clarify_val,
        r$diff_ratio)
  }
  add("")
}

# Benchmark section
add("BENCHMARKS  (n.sims=%d, %d reps, median elapsed)", N_BENCH, N_REPS)
add("--------------------------------------------------------------")
add("%-40s %8s %8s %6s", "Model", "post(s)", "clar(s)", "ratio")
add("--------------------------------------------------------------")
for (i in seq_len(nrow(bench_df))) {
  r <- bench_df[i, ]
  add("%-40s %8.3f %8.3f %6.2f", r$model, r$post_s, r$clarify_s, r$ratio)
}

add("")

# Summary section
add("SUMMARY")
add("--------------------------------------------------------------")
add("Total metric comparisons:   %d", test_count)
add("Passed:                     %d", pass_count)
add("Failed:                     %d", fail_count)
add("Errors (test did not run):  %d", error_count)
add("H4 (expected divergence):   scored separately above")
add("H5 (expected divergence):   scored separately above")
add("================================================================")

writeLines(lines, report_path)
cat(sprintf("Report written to %s\n", report_path))

# Print failures to console
if (fail_count > 0) {
  cat("\n--- FAILED COMPARISONS ---\n")
  fails <- all_df[!all_df$pass, ]
  print(fails, row.names = FALSE)
}

if (error_count > 0) {
  cat(sprintf("\nWARNING: %d test(s) encountered errors.\n", error_count))
}

cat(sprintf("\nDone: %d PASS / %d FAIL / %d ERROR\n",
            pass_count, fail_count, error_count))
