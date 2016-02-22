# Input y and z
y <- c(16/(58+16), 9/(9+90), 10/(10+48), 13/(13+57), 19/(19+103), 20/(20+57), 18/(18+86), 17/(17+112), 35/(35+273), 55/(55+64))
z <- c(12/(12+113), 1/(1+18), 2/(2+14), 4/(4+44), 9/(9+208), 7/(7+67), 9/(9+29), 8/(154+8))
mean <- seq(.0005, .5, .0005)
var <- seq(.00005, .05, .00005)

#Simulate posterior distribution of y and z
posterior_y <- matrix(rep(0, 1000000), 1000, 1000)
posterior_z <- matrix(rep(0, 1000000), 1000, 1000)
for (i in 1:1000){
    for (j in 1:floor(mean[i] * (1 - mean[i]) / 0.00025)){
        alpha <- mean[i] * mean[i] * (1 - mean[i]) / var[j] - mean[i]
        beta <- mean[i] * (1 - mean[i]) * (1 - mean[i]) / var[j] + mean[i] - 1
        posterior_y[i, j] <- 1 / (mean[i] * (1 - mean[i])) * prod(dbeta(y, alpha, beta))
        posterior_z[i, j] <- 1 / (mean[i] * (1 - mean[i])) * prod(dbeta(z, alpha, beta))
    }
}

posterior_y_m<-apply(posterior_y, 1, sum)
posterior_z_m<-apply(posterior_z, 1, sum)
alpha_y <- rep(0, 1000)
alpha_z <- rep(0, 1000)
beta_y <- rep(0, 1000)
beta_z <- rep(0, 1000)

for(x in 1:1000){
    i <- sample(1000, 1, prob = posterior_y_m)
    j <- sample(1000, 1, prob = posterior_y[i,])
    alpha_y[x] <- mean[i] * mean[i] * (1 - mean[i]) / var[j] - mean[i]
    beta_y[x] <- mean[i] * (1 - mean[i]) * (1 - mean[i]) / var[j] + mean[i] - 1
}

for(x in 1:1000){
    i <- sample(1000, 1, prob = posterior_z_m)
    j <- sample(1000, 1, prob = posterior_z[i,])
    alpha_z[x] <- mean[i] * mean[i] * (1 - mean[i]) / var[j] - mean[i]
    beta_z[x] <- mean[i] * (1 - mean[i]) * (1 - mean[i]) / var[j] + mean[i] - 1
}
plot(alpha_y, beta_y)
plot(alpha_z, beta_z)

miu_y <- alpha_y / (alpha_y + beta_y)
miu_z <- alpha_z / (alpha_z + beta_z)
diff <- miu_y - miu_z
hist(diff)
mean(diff)
