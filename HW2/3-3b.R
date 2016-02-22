#problem 3.3b
sample_c <- 1.013 + (0.24 / sqrt(32)) * rt(1000, 31)
sample_t <- 1.173 + (0.20 / sqrt(36)) * rt(1000, 35)
diff <- sample_t - sample_c
hist(diff, xlab="t - c", yaxt="n",breaks=seq(-0.1,0.4,0.02), cex=2)
quantile(diff, c(.025, .975))