library(triangle)

# The value of theta
thetas <- seq(0, 1, 0.001)

# Set prior P(theta), simulate the witch's hat distribution
prior_density <- 0.5 * dunif(thetas, 0, 1) + 0.5 * dtriangle(thetas, 0.385, 0.585)
names(prior_density) <- thetas

# Plot the Prior dist of theta
plot(prior_density ~ thetas, ylim = c(0, 6), type = 'l')

# The density of y|theta
likelihood <- lapply(thetas, function(x) {dbinom(437, size = 980, prob = x)})
likelihood <- unlist(likelihood)
names(likelihood) <- thetas

# The marginal dist of y
marginal_y <- sum(prior_density * likelihood)

# Posterior dist
post_theta_prob <- prior_density * likelihood / marginal_y
post_theta_density <- post_theta_prob / (0.001 * sum(post_theta_prob))
names(post_theta_prob) <- thetas
names(post_theta_density) <- thetas

# Plot of posterior density
plot(post_theta_density ~ thetas, ylim = c(0, 1.1 * max(post_theta_density)), type = 'l')


# posterior mean and variance
mean <- sum(thetas * post_theta_prob)
variance <- sum(thetas ^ 2 * post_theta_prob) - mean ^ 2

# posterior quantile function
quantile <- function(percentage) {
    cummu_prob_1 <- 0
    cummu_prob_2 <- 0
    for (theta in thetas) {
        cummu_prob_1 <- cummu_prob_1 + post_theta_prob[as.character(theta)]
        if ((cummu_prob_1 > percentage) && (cummu_prob_2) < percentage ) {
            print(theta)
            break
        } else {
            cummu_prob_2 <- cummu_prob_2 + post_theta_prob[as.character(theta)]
        }
    }
}

# Posterior median
quantile(0.5)

# 95% central posterior interval
quantile(0.025)
quantile(0.975)

