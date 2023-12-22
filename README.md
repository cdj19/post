# post
postestimation functions for regression in R

## description
The functions ‘post’ and ‘postSim’ are used to generate post-estimation quantities of interest from regression models using base R functions ‘lm’ and ‘glm’, as well as ‘svyglm’ from the ‘survey’ package and ‘polr’ from the ‘MASS’ package. ‘postSim’ is adapted from the `sim` function in the `arm` package, and takes draws from the sampling distribution of the estimated parameters. ‘post’ generates quantities of interest (e.g., first differences) and uses the simulated values from ‘postSim’ to calculate confidence bounds for these quantities. As a default, ‘post’ uses the so-called “observed value approach” in which the quantity is calculated for each observation and then averaged over the observations (as opposed to fixing independent variables at central tendencies or other researcher-chosen values). The user can fix values of subsets of model variables using the ‘holds’ option.

## variables
post(model,x1name=NULL,x1vals=NULL,x2name=NULL,x2vals=NULL,holds=NULL,n.sims=1000,cut=NULL,quantiles=c(.025,.975),did=NULL,weights=NULL, digits=2)

- model: an object of class ‘lm’, ‘glm’, ‘polr’, or ‘svyglm’.
- x1name: name of first variable, in quotes, that one wishes to vary for purposes of calculating predicted values or first differences.
- x1vals: a numeric vector containing a set of values of x1name at which to calculate predicted values of the dependent variable. A first difference is also calculated automatically using the first and last value in the vector.
- x2name: name of second variable, in quotes, that one wishes to vary for purposes of calculating predicted values or first differences.
- x2vals: a numeric vector containing a set of values of x2name at which to calculate predicted values of the dependent variable. One predicted value will be calculated for each combination of values in x1vals and x2vals: one conditional first difference for x1name will be calculated at each value in x2vals.
- holds: a list of variables to hold constant at specified values. (e.g., `holds = list(age = 35)`).
- n.sims: the number of draws to take from the sampling distribution.
- cut: an integer k. Used only for ‘polr’ models. When specified, ‘post’ returns the probability that Y exceeds its kth value at the specified values of x1name and x2name, as well as the associated first differences in this probability.
- quantiles: a numeric vector containing decimal values at which to calculate percentiles of the distribution of predicted values and first differences (confidence bounds) (e.g., for a 95% CI: `quantiles = c(0.025,0.975)`).
- did: a numeric vector containing two values from x2vals ordered from low to high. ‘post’ will calculate the difference in the first differences of x1name at these values of x2name with confidence bounds. When left unspecified, ‘post’ will use the minimum and maximum values of x2vals.
- weights: a numeric vector of weights to use in calculating the average over observations.
- digits: an integer specifying the number of decimal places to include in the output.

## output values
`est`: an array containing the requested quantities of interest (point estimates and first differences)
`did`: the difference in first differences (when relevant, otherwise empty)
`sims`: the simulated values of quantities of interest used to generate confidence bounds
`model`: the model class
`link`: the model link function
`quantiles`: the quantiles requested for confidence bounds
`call`: the user's call to `post`
