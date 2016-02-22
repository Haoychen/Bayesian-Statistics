library(mvtnorm)

# Input the data
x <- c(-0.86, -0.30, -0.05, 0.73)
n <- rep(5, 4)
y <- c(0, 1, 3, 5)
theta <- y / n

# Use logistic regression to calculate alpha and beta, nothing different with example in book
logis_regression <- glm(theta ~ x, family = binomial)
summary(logis_regression)

# Compute the posterior distribution
miu<-c(0,10)
alpha <- seq(-5, 10, 0.05)
beta <- seq(-10, 40, 0.05)
sigma <- matrix(c(4, 10, 10, 100), nrow = 2)
posterior <- matrix(0,length(alpha),length(beta))

for (i in 1:length(alpha)){
    for (j in 1:length(beta)){
        posterior[i,j]<-dmvnorm(c(alpha[i],beta[j]),miu,sigma)*prod(dbinom(y,n,exp(alpha[i]+beta[j]*x)/(1+exp(alpha[i]+beta[j]*x))))
    }
}

# Draw contour plot
contour(alpha, beta, posterior, xlim = c(-2,4), ylim = c(-1,25), xlab = "alpha", ylab ="beta")

# Draw scatterplot
posterior <- posterior/sum(posterior)	
randalpha <- rep(0,1000)
randbeta <- rep(0,1000)
alpha_y <- rowSums(posterior)
alpha_ycum <- cumsum(alpha_y)

for(i in 1: 1000){
    ind <- which(alpha_ycum>=runif(1,0,1))[1]
    randalpha[i] <- alpha[ind]
    beta_alpha_y <- posterior[ind,]/alpha_y[ind]
    beta_alpha_ycum <- cumsum(beta_alpha_y)
    randbeta[i] <- beta[which(beta_alpha_ycum>=runif(1,0,1))[1]]
}

plot(randalpha,randbeta, xlim = c(-2,6), ylim = c(0,30), xlab = "alpha", ylab ="beta")


# LD50
LD50 <- randalpha / randbeta
hist(LD50)