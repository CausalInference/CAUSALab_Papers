############################
## Example of quadratic forms
## Author: Elizabeth Diemer
###########################

## needed function
inv.logit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

# generate possible a values
# discrete for now for laziness
a <- seq(0,7, by = .1)
beta_y23 <- seq(1, 7.9, by = 0.5)
values <- list(a, beta_y23)
eqn_values <- expand.grid(values, stringsAsFactors = FALSE)
names(eqn_values) <- c("a", "beta_y23")
eqn_values <- as.data.frame(eqn_values)
## generate pred prob of y from a and beta_y23
eqn_values$predy <- inv.logit(-1.31 - 0.2675*eqn_values$a +log(eqn_values$beta_y23)*eqn_values$a^2)

##graph predicted probabilities
eqn_values$beta_y23 <- as.character(eqn_values$beta_y23)
ggplot(data = eqn_values, 
       aes(x = a, y = predy, color = beta_y23)) +
  geom_line(aes(linetype = beta_y23)) + labs(x = "A", 
                                             y = "Predicted Probability of Y", 
                                             color = "Quadratic Term",
                                             linetype = "Quadratic Term") 
