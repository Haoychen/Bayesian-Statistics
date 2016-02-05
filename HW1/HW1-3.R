library(triangle)

# The value of theta
theta <- seq(0, 1, 0.001)

# Set prior P(theta), simulate the witch's hat distribution, 40%
prior_density <- 0.5 * dunif(theta, 0, 1) + 0.5 * dtriangle(theta, 0.385, 0.585)
names(prior_density) <- theta

# The density of y|theta
likelihood <- lapply(theta, function(x) {dbinom(437, size = 980, prob = x)})
likelihood <- unlist(likelihood)
names(likelihood) <- theta

# The marginal dist of y
marginal_y <- sum(prior_density * likelihood)

# Posterior dist
post_theta_density <- prior_density * likelihood / marginal_y

# Plot of posterior density
plot(post_theta_density ~ theta)

