# post
Post-estimation functions for regression in `R`

## Description

The functions `post` and `postSim` generate post-estimation quantities of interest from regression models fit with base R functions `lm` and `glm`, as well as `svyglm` from the `survey` package and `polr` from the `MASS` package.

`post` computes predicted values, first differences, and differences-in-differences on the response scale, with confidence intervals. It uses the "observed value" approach by default: the quantity of interest is calculated for every observation in the estimation sample and then averaged, rather than fixing covariates at means or other researcher-chosen values. Subsets of variables can be held constant at specific values using the `holds` option.

Two inference methods are available:

- **Simulation** (default): `postSim` draws from the sampling distribution of the estimated parameters and `post` uses the empirical distribution of simulated quantities to construct confidence intervals.
- **Delta method**: `post` analytically computes standard errors via the Delta method and constructs normal-approximation confidence intervals, without simulation.

## Supported models and link functions

| Model class | Package | Supported link functions |
|---|---|---|
| `lm` | base R | identity |
| `glm` | base R | identity, logit, probit, log, cloglog |
| `svyglm` | `survey` | identity, logit, probit, log, cloglog |
| `polr` | `MASS` | logistic, probit, cloglog |

## `post`

```r
post(model, x1name = NULL, x1vals = NULL, x2name = NULL, x2vals = NULL,
     holds = NULL, n.sims = 1000, cut = NULL, quantiles = c(.025, .975),
     did = NULL, weights = NULL, digits = 2,
     dist = c("normal", "t"), method = c("simulation", "delta"))
```

### Parameters

- **`model`**: a fitted model object of class `lm`, `glm`, `svyglm`, or `polr`.
- **`x1name`**: character string naming a variable in the model to vary. Predicted values are computed at each value in `x1vals`, and a first difference is computed between the first and last values.
- **`x1vals`**: a vector of values at which to evaluate `x1name`. For factor variables, supply character values matching the factor levels (e.g., `x1vals = c("low", "high")`).
- **`x2name`**: character string naming a second variable to vary. When specified along with `x1name`, predicted values are computed for every combination of `x1vals` and `x2vals`, a conditional first difference for `x1name` is computed at each value in `x2vals`, and a difference-in-differences is computed.
- **`x2vals`**: a vector of values at which to evaluate `x2name`. Factor levels are supported.
- **`holds`**: a named list of variables to hold constant at specified values (e.g., `holds = list(age = 35, region = "west")`). All other variables remain at their observed values.
- **`n.sims`**: integer; the number of draws from the sampling distribution. Used only when `method = "simulation"`.
- **`cut`**: integer *k*. For `polr` models only. When specified, `post` returns P(Y > *k*th category) instead of probabilities for each category, along with associated first differences.
- **`quantiles`**: numeric vector of quantiles for confidence bounds (e.g., `c(0.025, 0.975)` for a 95% CI).
- **`did`**: a vector of two values from `x2vals` specifying which pair to use for the difference-in-differences calculation. Defaults to the first and last values of `x2vals`.
- **`weights`**: numeric vector of observation weights for the averaging step. Must be the same length as the estimation sample.
- **`digits`**: integer; number of decimal places in the output.
- **`dist`**: `"normal"` (default) or `"t"`. Controls the sampling distribution used by `postSim`. When `"t"`, draws are taken from a multivariate *t* distribution with residual degrees of freedom instead of the multivariate normal. Used only when `method = "simulation"`.
- **`method`**: `"simulation"` (default) or `"delta"`. When `"delta"`, standard errors are computed analytically via the Delta method and confidence intervals use the normal approximation. No simulation is performed, so `n.sims` and `dist` are ignored.

### Output

`post` returns an S4 object with the following slots:

- **`est`**: array of point estimates and confidence bounds. Dimensions depend on the call:
  - No `x1name`: a single predicted value with CI.
  - `x1name` only: one row per value of `x1vals`, plus a first difference row.
  - `x1name` + `x2name`: a 3D array (x1 values + first difference) x (quantiles) x (x2 values).
  - For `polr` without `cut`: an additional dimension for each outcome category.
- **`did`**: array containing the difference-in-differences estimate with CI, or `NULL` when `x2name` is not specified.
- **`sims`**: matrix of simulated quantities of interest (when `method = "simulation"`), or `NULL` (when `method = "delta"`).
- **`model`**: character string with the model class.
- **`link`**: the link function used.
- **`quantiles`**: the quantiles used for confidence bounds.
- **`call`**: the original call to `post`.

## `postSim`

```r
postSim(object, n.sims = 1000, dist = c("normal", "t"))
```

`postSim` draws from the sampling distribution of the model's estimated parameters. It is adapted from the `sim` function in the `arm` package.

### Parameters

- **`object`**: a fitted model object of class `lm`, `glm`, `svyglm`, or `polr`.
- **`n.sims`**: integer; the number of draws.
- **`dist`**: `"normal"` (default) or `"t"`. When `"t"`, draws are taken from a multivariate *t* distribution with residual degrees of freedom.

### Output

- For `lm`, `glm`, and `svyglm` models: a `postSim` object with slots `coef` (matrix of coefficient draws) and `sigma` (vector of residual SD draws).
- For `polr` models: a `postSim.polr` object with slots `coef` (matrix of coefficient draws) and `zeta` (matrix of cutpoint draws).

## Examples

```r
source("postSim.R")
source("post.R")

# Linear model — predicted values and first difference for wt
fit <- lm(mpg ~ wt + hp, data = mtcars)
post(fit, x1name = "wt", x1vals = c(2, 3, 4))

# GLM logit — simulation-based (default)
fit <- glm(vs ~ wt + hp, data = mtcars, family = binomial())
post(fit, x1name = "wt", x1vals = c(2, 4), n.sims = 5000)

# GLM logit — Delta method
post(fit, x1name = "wt", x1vals = c(2, 4), method = "delta")

# Difference-in-differences
post(fit, x1name = "wt", x1vals = c(2, 4),
     x2name = "hp", x2vals = c(100, 200))

# Holding a variable constant
post(fit, x1name = "wt", x1vals = c(2, 4), holds = list(hp = 150))

# Ordered logit — category probabilities
library(MASS)
mtcars$gear_o <- ordered(mtcars$gear)
fit <- polr(gear_o ~ wt + hp, data = mtcars, method = "logistic")
post(fit, x1name = "wt", x1vals = c(2, 4))

# Ordered logit — P(Y > 3) using cut
post(fit, x1name = "wt", x1vals = c(2, 4), cut = 1)

# Factor variables
mtcars$am_f <- factor(ifelse(mtcars$am == 1, "manual", "auto"))
fit <- glm(vs ~ am_f + wt, data = mtcars, family = binomial())
post(fit, x1name = "am_f", x1vals = c("auto", "manual"))

# Custom weights and quantiles
post(fit, x1name = "wt", x1vals = c(2, 4),
     weights = runif(nrow(mtcars)),
     quantiles = c(0.05, 0.95))

# Multivariate t draws
post(fit, x1name = "wt", x1vals = c(2, 4), dist = "t")
```

## Dependencies

- `MASS` (for `mvrnorm` and `polr`)
- `survey` (optional, for `svyglm` models)
