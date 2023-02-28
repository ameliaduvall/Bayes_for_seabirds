## Introduction to Bayesian Analysis for Seabird Research
## Workshop at PSG 2023 Annual Meeting


######### Libraries and data ##############################
# Libraries (install packages if needed)
library(rstan)
library(here)

# Define the murre dataset and create a plot
# Data is illustrative only
murre.df <- data.frame(
  year=c(1990,1993,1995,1997,
         2002,2002,2002,
         2003,2003,2003,2004,2004,2005,2005,2005,2006,2007,2007,2008,
         2009,2009,2010,2010,2011,2011,2012,2012,2012,2013,2013,
         2014,2014,2015,2016,2016,2017,2017,2018,2018,2019,2020,2020,2020),
  count=c(14540,10850,14710,17760,
          10960,13180,12110,11080,10680,12530,12340,11385,
          8690,9450,9720,8500,
          8710,8690,11130,9430,9700,14190,9850,8750,8890,
          10090,12640,10960,10940,12520,10740,9830,9960,
          6470,7140,7860,6890,7720,8910,7990,8500,9010,8020)
)

# Plot
plot(count ~ year, data=murre.df, type="p", pch=21, cex=2,lwd=2,
     ylab="Count", xlab="Year", bg="white")

###########################################################

######### Prior definition examination ####################
# Using a poisson GLM, need intercept and slope
# In the model we will re-scale the year variable
# so that t = 0 corresponds to year 2005 (mid-point of data)

# Intercept (B0)
# Murre colonies range from 10's to millions of birds
# Define our prior to 
# - approach 0 closer to B0 = log(1) = 0
# - approach 0 closer to B0 = log(1e7) = 16.12
B0.prior <- function(x){dnorm(x, mean=8.0, sd=4)}

# Plot the probability distributions
curve(B0.prior,from = -2.5, to = 20, type="l",
     ylab="Probability density", 
     xlab=expression(paste("Parameter: ", beta[0])), lwd=2)
abline(h=0, col=rgb(0,0,0,0.2))
text(x=15, y=0.08,
     labels=expression(paste("Prior: ", beta[0], " ~ Normal(8,4)")))
axis(side=3, at=log(c(0.1,1,10,100,1e3,1e4,1e5,1e6,1e7)),
     labels=c(0.1,1,10,100,1000,1e4,1e5,1e6,1e7),
     line=-0)

# Slope (B1)
# A consistent 10-fold increase/decrease seem unlikely
# without any other knowledge, might expect constant
# centre prior on zero with 10-fold change @ 3SD
B1.prior <- function(x){dnorm(x, mean=0, sd=0.8)}

# Plot the probability distributions
curve(B1.prior,from = -3, to = 3, type="l",
      ylab="Probability density", 
      xlab=expression(paste("Parameter: ", beta[1])), lwd=2)
abline(h=0, col=rgb(0,0,0,0.2))
text(x=2, y=0.4,
     labels=expression(paste("Prior: ", beta[1], " ~ Normal(0,0.8)")))
axis(side=3, at=log(c(0.1,0.2,0.5,1,2,5,10)),
     labels=c(0.1,0.2,0.5,1,2,5,10),line=-0,)

###########################################################

############## STAN GLM USING RSTANARM ####################
library(rstanarm)
library(shinystan)

# Offset year so that 2005 is t=0
murre.df$year.off <- murre.df$year - 2005

# The simple way
murre_stanarm <- stan_glm(count ~ year.off, data = murre.df,
                       family = poisson(link="log"),
                       prior_intercept = normal(8,4),
                       prior = normal(0,0.8),
                       iter = 2000,
                       warmup = 1000,
                       chains = 4)

# examine coefficients
print(murre_stanarm, digits=4)

# create 95% credible interval
posterior_interval(murre_stanarm, prob=0.95)

# view diagnostic plots
launch_shinystan(murre_stanarm)

############################################################

############## STAN MODEL FITTING USING BASE STAN #########
# see poisson_glm_psg.stan

####
# format the data as inputs to STAN in a list
mod.in = list(N_obs=nrow(murre.df), 
              Year=murre.df$year.off, 
              Murres=murre.df$count)
####
# Fit model
murre.model <- stan(
  file = here("code", "poisson_glm_psg.stan"), # Model definition
  data = mod.in,               # List of data inputs
  chains = 4,                  # Number of MCMC chains
  warmup = 4000,               # Number of warm-up iterations
  iter = 8000,                 # Total number of iterations
  verbose = FALSE,             # Verbose?
  refresh=200)                 # How often to report progress
####

#####################################################################

################### MODEL CHECKS AND EVALUATION #####################
# There are a set of criteria that we want to look out for in
# terms of whether the model converged on a solution
# The statistic R-hat is a measure of chain mixing, and we ideally
# want that to be < 1.02, and closer to zero the better
summary(murre.model)$summary

print(murre.model, pars=c("beta0","beta1"), probs=c(.1,.5,.9))

# Examine posterior distributions
mod.params <- extract(murre.model, pars=c("beta0","beta1"))

# Plot histograms
par(mfrow=c(1,2))
hist(mod.params$beta0, breaks=100, main="Intercept", xlab="",
     freq=F)
hist(mod.params$beta1, breaks=100, main="Slope", xlab="",
     freq=F)

# Or view things using the nice shinystan
launch_shinystan(murre.model)

#################################################################

################### POSTERIOR PREDICTION ########################
# In the model code, I included a prediction for the expected count
# in 1998 in the generated quantities block.
# We can extract this in the same way as for model parameters
c1998 <- extract(murre.model, pars=c("c1998"))

hist(c1998$c1998, breaks=100,
     xlab="Predicted abundance in 1998",
     ylab="Probability density", freq=F,
     main="")

# NOTE: the poisson distribution is very rarely applicable
# to many counts, it was simply used here as an example
# within the stan code file there are additional lines 
# for exploring the negative binomial version which you
# are welcome to explore

################################################################

