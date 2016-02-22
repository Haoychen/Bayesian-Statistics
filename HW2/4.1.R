#problem 4.1 b)
library(stats)
y <- c(43,44,45,46.5,47.5)
first_derivative <- function(theata){
    2 * sum((y - theata) / (1 + (y - theata) ^ 2))
}
mode <- uniroot(first_derivative, lower = -100, upper = 100)[1]
mode <- as.numeric(mode)

# 4.1c


#plot the exact density
theta <- seq(0,100,0.001)
y <- c(43,44,45,46.5,47.5)
posterior <-rep(0, length(theta))
for(i in 1:length(theta)){
    posterior[i] <- prod(dcauchy(y, location = theta[i]))
}
plot(theta, posterior / (sum(posterior) * 0.001), xlim = c(40,50), ylim = c(0,0.5), lty = 2, type = 'l')
curve(dnorm(x, mean = mode, sd = sqrt(1 / 1.37489)), add = T)

