# b
y <- c(24,25,31,31,22,21,26,20,16,22)
alpha <- seq(mean(y) - 4 * sd(y), mean(y) + 4*sd(y), length.out = 1000)
beta <- seq(-4, 4, length.out = 1000)
prior<-matrix(rep(0,1000000), 1000, 1000)
for (i in 1:1000){
    for (j in 1:1000){
        if (alpha[i] + 12 * beta[j] > 0){
            prior[i, j] <- exp(-(alpha[i] - mean(y)) ^ 2 * 0.5 / var(y) - beta[j] ^ 2 / 2)
        }
    }
}
contour(alpha, beta, prior, xlim = c(10, 40), ylim = c(-4, 4), xlab = "alpha", ylab = "beta", main = "contour of informative prior distribution")

# e
t <- seq(1:10)
fit <- lm(y ~ t)
summary(fit)


# f
posterior <- matrix(rep(0, 1000000), 1000, 1000)
for (i in 1:1000){
    for (j in 1:1000){
        if (alpha[i] + 12 * beta[j] > 0){
            # simulate posterior
            posterior[i, j] <- prod(exp(-(alpha[i] + beta[j] * t)) * (alpha[i] + beta[j] * t) ^ y)
        }
    }
}
contour(alpha, beta, posterior, xlim=c(10,40), ylim=c(-4,4), xlab="alpha", ylab="beta", main="contour of joint posterior density")

#take 1000 draws from the joint posterior density
post_alpha <- apply(posterior, 1, sum)
alpha_sim<-rep(0, 1000)
beta_sim<-rep(0, 1000)
for (k in 1:1000){
    i <- sample(1000, 1, prob = post_alpha)
    j <- sample(1000, 1, prob=posterior[i,])
    alpha_sim[k] <- alpha[i]
    beta_sim[k] <- beta[j]
}
plot(alpha_sim, beta_sim, xlim=c(10,40), ylim=c(-4,4), xlab="alpha", ylab="beta", main="simulation from the joint posterior density")

# g
hist(alpha_sim + 11 * beta_sim, breaks = 20, xlab="expected number", main="posterior density for the expected number of fatal accidents")

# h
sim <- rpois(1000, alpha_sim + 11 * beta_sim)
quantile(sim, c(0.025,0.975))

# i
contour(alpha, beta, prior, xlim=c(10,40), ylim=c(-4,4),xlab="alpha",ylab="beta")
points(alpha_sim, beta_sim)

# a
theta <- -seq(-10, 10, 0.01)
p <- pnorm(theta/4 + 0.675 * sqrt(5/4)) - pnorm(theta / 4 - 0.675 * sqrt(5/4))
p_theta <- 1 / sqrt(8 * 3.14159265) * exp(-theta ^ 2 /8)
sum(0.01 * p * p_theta)

# c
plot(theta, p, type = 'l')


# b
y <- c(16,9,10,13,19,20,18,17,35,55)
n <- c(58,90,48,57,103,57,86,112,273,64) + y
r <- y/n
alpha <- (mean(r) ^ 2 * (1 - mean(r))) / var(r) - mean(r)
beta <- (mean(r) * (1 - mean(r)) ^ 2) / var(r) + mean(r) - 1

a <- seq(-2.0, -0.5, length.out=1000)
b <- seq(1.0, 4.5, length.out=1000)
prior <- matrix(0, 1000, 1000)
posterior <- matrix(0, 1000, 1000)
for (i in 1:1000){
    for (j in 1:1000){
        alpha_sim <- exp(a[i] + b[j]) / (exp(a[i]) + 1)
        beta_sim <- exp(b[j]) / (exp(a[i]) + 1)
        prior[i, j] <- (alpha_sim + beta_sim) ^ (-2.5) * alpha_sim * beta_sim
        posterior[i, j] <- prior[i, j] * prod(beta(alpha_sim + y, beta_sim + n - y) / beta(alpha_sim, beta_sim))
    }
}
contour(a, b, posterior/max(posterior), xlab="log(alpha/beta)", ylab="log(alpha+beta)", main="marginal posterior density")

post_alpha <- apply(posterior, 1, sum)
alpha_sim = rep(0,1000)
beta_sim = rep(0,1000)
for(i in 1:1000){
    loc <- sample(1000, 1, prob=post_alpha)
    a_sim <- a[loc]
    post_beta <- posterior[loc,]
    b_sim <- sample(b, 1, prob=post_beta)
    alpha_sim[i] <- exp(a_sim + b_sim) / (exp(a_sim) + 1)
    beta_sim[i] <- exp(b_sim) / (exp(a_sim) + 1)
}
plot(alpha_sim, beta_sim, main="simulations from the joint posterior distribution")


# c
library(plotrix)

theta <- matrix(rep(0, 10 * 1000), 10, 1000)
for(i in 1:10)
    for(j in 1:1000){
        {
            theta[i,j]=rbeta(1, alpha_sim[j] + y[i], beta_sim[j] + n[i] - y[i])
        } 
    }
mean <- rowMeans(theta)
median <- apply(theta, 1, median)
lower <- vector()
upper <- vector()

for(i in 1:10){
    lower[i] <- quantile(theta[i,],0.025)
    upper[i] <- quantile(theta[i,],0.975)
}

theta_raw <- y / n
plotCI(theta_raw, median, ui=upper, li=lower, xlab="observed rate", ylab="95% condidence interval of simulation")
abline(a=0,b=1)

# d
mean(lower)
mean(upper)

# e
y_new<-seq(0, 100, 0.1)
theta_avg<-colMeans(theta)
post_pred<-rep(0, 101)
for (i in 1:1000){
    post.pred<-post_pred + dbinom(y_new, 100, theta_avg[i]) / 1000
}

prob <- 0
CI <- c(0, 0)
for (i in y_new){
    prob <- prob + post_pred[i + 1]
    if(prob <= 0.025) CI[1] <- i
    if(prob >= 0.975) {
        CI[2] <- i
        break
    }
}
CI


