library(foreign)
library(plyr)
library(dplyr)
library(openintro)

# read data and select useful features
elect.data <- read.dta('pew_research_center_june_elect_wknd_data.dta')
Obama_votes <- read.csv('2008ElectionResult.csv')
elect.data <- elect.data[,c('ideo', 'state')]
Obama_votes <- Obama_votes[, c('state', 'vote_Obama_pct')]

# Exclude Alaska, Hawaii and Washington DC and rows with NA
elect.data <- filter(elect.data, !(state %in% c('alaska', 'washington dc', 'hawaii')))
elect.data <- elect.data[complete.cases(elect.data),]

#Transfer state name to state abbreviation
elect.data$state <- state2abbr(elect.data$state)
Obama_votes$state <- state2abbr(Obama_votes$state)

# split the election data by state, calculate the proportion of very liberal people 
# and combine together

liberal_proportion <- ddply(elect.data, 'state', function(dataframe){
    liberal_people <- dataframe[dataframe$ideo == 'very liberal',]
    proportion <- nrow(liberal_people) / nrow(dataframe)
    data.frame(liberal_people = nrow(liberal_people), total_people = nrow(dataframe), state_proportion = proportion)
})

# Combine liberal_probortion dataset and Obama_votes dataset by state abbreviation
liberal_proportion_vs_Obama_vote <- merge(liberal_proportion, Obama_votes, by = 'state')

# Plot the liberal proportion vs Obama vote
plot(liberal_proportion_vs_Obama_vote$vote_Obama_pct, liberal_proportion_vs_Obama_vote$state_proportion,
     ylab = 'Liberal Proportion', xlab = 'Obama Vote Share', type = 'n')
text(liberal_proportion_vs_Obama_vote$vote_Obama_pct, liberal_proportion_vs_Obama_vote$state_proportion,
     labels = liberal_proportion_vs_Obama_vote$state, cex = 0.7, pos = 3)

# Estimate the prior
prior_mean <- mean(liberal_proportion$state_proportion)
prior_var <- var(liberal_proportion$state_proportion)

alpha <- (prior_mean ^ 2) * (1 - prior_mean) / prior_var - prior_mean
beta <- alpha * (1 - prior_mean) / prior_mean

# Posterior Mean for each state
liberal_proportion_vs_Obama_vote$posterior_mean <- (alpha + liberal_proportion$liberal_people) / (alpha + beta + liberal_proportion$total_people)

# Plot the posterior mean vs Obama vote
plot(liberal_proportion_vs_Obama_vote$vote_Obama_pct, liberal_proportion_vs_Obama_vote$posterior_mean,
     ylab = 'Posterior Mean', xlab = 'Obama Vote Share', type = 'n')
text(liberal_proportion_vs_Obama_vote$vote_Obama_pct, liberal_proportion_vs_Obama_vote$posterior_mean,
     labels = liberal_proportion_vs_Obama_vote$state, cex = 0.7, pos = 3)

# Plot the liberal Mean vs number of respondents
plot(liberal_proportion_vs_Obama_vote$total_people, liberal_proportion_vs_Obama_vote$state_proportion,
     ylab = 'Liberal Proportion', xlab = 'Number of Respondent', type = 'n')
text(liberal_proportion_vs_Obama_vote$total_people, liberal_proportion_vs_Obama_vote$state_proportion,
     labels = liberal_proportion_vs_Obama_vote$state, cex = 0.7, pos = 3)


# Plot the liberal Mean vs Number of Respondents
plot(liberal_proportion_vs_Obama_vote$total_people, liberal_proportion_vs_Obama_vote$posterior_mean,
     ylab = 'Poesterior Mean', xlab = 'Number of Respondent', type = 'n')
text(liberal_proportion_vs_Obama_vote$total_people, liberal_proportion_vs_Obama_vote$posterior_mean,
     labels = liberal_proportion_vs_Obama_vote$state, cex = 0.7, pos = 3)
