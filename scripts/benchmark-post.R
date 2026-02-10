# Minimal benchmark to compare post() execution paths after bug fixes.
# Run with: Rscript scripts/benchmark-post.R

source("postSim.R")
source("post.R")

fit <- glm(vs ~ wt + hp + factor(cyl), data = mtcars, family = binomial())

set.seed(123)
base_time <- system.time(
  replicate(50, post(fit, n.sims = 100, digits = 6), simplify = FALSE)
)

set.seed(123)
hold_time <- system.time(
  replicate(50, post(fit, holds = list(cyl = 6), n.sims = 100, digits = 6), simplify = FALSE)
)

print(base_time)
print(hold_time)
