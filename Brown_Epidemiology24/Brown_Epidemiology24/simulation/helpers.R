expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit <- function(x) {
  return(log(x/(1-x)))
}

sandwich1 <- function(object, ...) {sandwich(object) * nobs(object) / (nobs(object) - 1)}
