#5.14
#b.
y<-c(16,9,10,13,19,20,18,17,35,55)
y<-c(58,90,48,57,103,57,86,112,273,64) + y

log_alpha_over_beta <- seq(3, 7, length.out=100) 
log_beta<-seq(-8, -1, length.out=100)
beta<-exp(log_beta)
alpha = exp(log_alpha_over_beta) * beta
I<-length(alpha)
J<-length(beta)
posterior<-matrix(rep(NA, I*J), I, J)

fi = 0.001
for(i in 1:I){
  for(j in 1:J){
    fact <- (alpha[i] * beta[j]) ^ (fi - 1) * exp(-fi * (alpha[i] + beta[j]))
    posterior[i, j] <- fact * prod(beta[j] ^ alpha[i] / (beta(alpha[i], y) * (beta[j] + 1) ^ (alpha[i] + y)))
  }
}
contour(log_alpha_over_beta,log_beta, posterior / max(posterior), drawlabels=F, xlab="log(alpha/beta)", ylab="log(beta)")

post.alpha <- apply(posterior, 1, sum)
alpha.sim <- vector()
beta.sim <- vector()

for(i in 1:1000)
{
  foo <- sample(length(alpha), 1, prob = post.alpha)
  alpha.sim[i] <- alpha[foo]
  post.beta <- posterior[foo,]
  beta.sim[i] <- sample(beta, 1, prob = post.beta)
}

plot(alpha.sim, beta.sim, xlim = c(0, 8), ylim = c(0, 0.1), xlab = "alpha", ylab = "beta")
plot(log(alpha.sim / beta.sim), log(beta.sim), xlim=c(1,12), ylim = c(-8, -1), xlab = "log(alpha/beta)", ylab = "log(beta)")

#e
y <- c(16,9,10,13,19,20,18,17,35,55)
n <- c(58,90,48,57,103,57,86,112,273,64) + y
theta <- matrix(rep(NA,10*1000),10,1000)
for(i in 1:10)
  for(j in 1:1000)
  {
    {
      theta[i,j] <- rgamma(1,alpha.sim[j] + y[i], beta.sim[j] + n[i])
    } 
  }

mean <- rowMeans(theta)
median <- apply(theta, 1, median)
medlower <- vector()
upper <- vector()

for(i in 1:10){
  lower[i] <- quantile(theta[i,],0.025)
  upper[i] <- quantile(theta[i,],0.975)
}

theta.raw <- y/n
library(plotrix)
plotCI(theta.raw, median, ui=upper, li=lower, xlab="observed rate", ylab="95% condidence interval of simulation")
abline(a=0, b=1)


#5.15

#a.
x <- read.table("meta.asc.txt", head = T)
x <- x[,-1]
y <- -log(x[,1]/(x[,2]-x[,1]))+log(x[,3]/(x[,4]-x[,3]))
sigma <- sqrt(1/x[,1]+1/(x[,2]-x[,1])+1/x[,3]+1/(x[,4]-x[,3]))

tau <- seq(0.001,2.0,0.001)
posterior <- vector()
mu <- vector()
for(i in 1:length(tau))
{
  mu.est <- sum(y / (sigma^2 + tau[i]^2)) / sum(1 / (sigma^2 + tau[i]^2))
  v.mu <- 1 / sum(1 / (sigma^2 + tau[i]^2))
  prior <- 1
  mu[i] <- rnorm(1, mean = mu.est, sd = sqrt(v.mu))
  posterior[i] <- prior * sqrt(v.mu) * prod((sigma^2 + tau[i]^2) ^ (-0.5) * exp(-(y - mu.est) ^ 2 / (2 * (sigma^2 + tau[i]^2))))
}

plot(tau, posterior, type='l', xlim=c(0, 0.6))

#b
theta <- matrix(rep(NA, length(y) * length(tau)), length(y), length(tau))
theta.mean <- matrix(rep(NA, length(y) * length(tau)), length(y), length(tau))
theta.sd <- matrix(rep(NA, length(y) * length(tau)), length(y), length(tau))
for(i in 1:length(y))
  for(j in 1:length(tau))
  {
    {
      mu <- sum(y / (sigma^2 + tau[j]^2)) / sum(1 / (sigma^2 + tau[j]^2))
      temp_mean <- (y[i] / sigma[i]^2 + mu / tau[j]^2) / (1 / sigma[i]^2 + 1 / tau[j]^2)
      v<-1/(1/sigma[i]^2+1/tau[j]^2)
      theta[i,j] <- rnorm(1,mean=temp_mean,sd=sqrt(v))
      theta.mean[i,j] <- temp_mean
      theta.sd[i,j] <- sqrt(v)
    } 
  }

plot(tau,theta.mean[1,], type='l', ylim=c(-1,0.5),ylab="Estimated Treatment Effect")
for(i in 1:21)
{
  lines(tau,theta.mean[(i+1),], col=i)
}

plot(tau,theta.sd[1,],type='l', ylim=c(0,1),ylab="Posterior Standard Deviation")
for(i in 1:21)
{
  lines(tau,theta.sd[(i+1),],col=i)
}


#c
median<-apply(theta,1,median)

plot(median,y,pch=20)
text(median,y,1:J,pos=2,cex=0.6)
abline(a=0,b=1)
abline(v=mean(theta))

#d
theta.new<-vector()
for(i in 1:length(tau)){
  foo<-sample(length(posterior),1,prob=posterior)
  tau.new<-tau[foo]
  mu.new<-sum(y/(sigma^2+tau.new^2))/sum(1/(sigma^2+tau.new^2))
  theta.new[i]<-rnorm(1,mean=mu.new,sd=tau.new)
}
hist(theta.new,main="simulation from posterior distributionof a new treatment effect",breaks=20)

#e
p0<-x$V2/x$V3
p1<-x$V4/x$V5
y.new<-vector()
for(i in 1:length(tau))
{
  index<-sample(c(0,1),1)
  loc<-sample(length(p0),1)
  if(index==0){ ps0<-p0[loc]
  ps1<-p0[loc]*exp(theta.new[i])
  }
  else {ps1<-p1[loc]
  ps0<-ps1/exp(theta.new[i])}
  y.new[i]<-rnorm(1,mean=theta.new[i],sd=sqrt((1/ps0+1/ps1+1/(1-ps0)+1/(1-ps1))/100))
}
hist(y.new,main="Histogram of simulation of y.new|y")


#6.6 b
test<-vector()
for(i in 1:1000){
  theta <- rbeta (1, 8, 14)
  y.rep <- rbinom (1, 1, theta)
  while(sum(y.rep==0) < 13){
    y.rep <- c(y.rep, rbinom(1, 1, theta))
  }
  n.rep <- length(y.rep)
  test[i] <- sum(y.rep[2:n.rep] != y.rep[1:(n.rep-1)])
}
hist (test, xlab="T (y-rep)", yaxt="n",
      breaks=seq(-.5, max(test) + .5), cex = 2)
#6.7
test <- vector()
for (i in 1:1000){
  theta <- rnorm (1,5.1,0.1)
  y.rep <- rnorm (100,theta,1)
  test[i] <- max(abs(y.rep))
}
hist(test, breaks = 20, xlab = "T(y_rep)")
lines(rep(8.1,2), c(0,1000))
mean(test > 8.1)

