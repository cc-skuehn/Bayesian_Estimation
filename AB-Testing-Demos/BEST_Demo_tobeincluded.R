# The complete R Source code of BEST Demo Notebook (all chunks combined)

#############
### Setup ###
#############
## Set path
setwd('~/git/analytics_scripts/cross_bu/AB-Testing-Demos/')

## Helper file for better plotting
source('~/git/analytics_scripts/cross_bu/AB-Testing-Demos/HelperFiles/Multiplot_fun.R')

## Libraries
library(ggplot2)

## Set random seed for reproducibility
set.seed(7)
#set.seed(70)

#######################
### Example problem ###
#######################
## Create samples
sample.1 <- rnorm(80, 100, 3)
sample.2 <- rnorm(100, 103, 7)

## Pool data for parameter estimation in the prior
## Reasoning: Null hypothesis is that the data/samples are drawn from the same distribution, that's why we can pool the samples under this hypothesis
pooled <- c(sample.1, sample.2)
## Inspect sample distributions
## Blue: Mean of the underlying "generating" distribution
## Darkgreen: Mean of the concrete sample
hs1 <- ggplot()+geom_histogram(aes(sample.1),binwidth = 1,fill='black')+xlim(min(pooled)-1,max(pooled)+1)+geom_vline(xintercept = 100,color='blue')+geom_vline(xintercept = mean(sample.1),color='darkgreen')+theme_bw()
hs2 <- ggplot()+geom_histogram(aes(sample.2),binwidth = 1,fill='black')+xlim(min(pooled)-1,max(pooled)+1)+geom_vline(xintercept = 103,color='blue')+geom_vline(xintercept = mean(sample.2),color='darkgreen')+theme_bw()
multiplot(hs1,hs2)

## Likelihood and prior
## Technically, the "evidence" is missing in the formulas, in Bayesian terminology: posterior = likelihood * prior / evidence
## But as evidence is constant and practically not computable it is omitted here
## Still, the scaling effect of the evidence is important as for larger sample sizes the prod() will be 0 or NaN (due to rounding)
## Use the BEST package in order to avoid this: install.packages("BEST"), this installs rjags as well (required for the MCMC sampling part)
## Likelihood
likelihood <- function(parameters){
  mu1=parameters[1]; sig1=parameters[2]; mu2=parameters[3]; sig2=parameters[4]
  prod(dnorm(sample.1, mu1, sig1)) * prod(dnorm(sample.2, mu2, sig2))
}
## Prior
prior <- function(parameters){
  mu1=parameters[1]; sig1=parameters[2]; mu2=parameters[3]; sig2=parameters[4]
  dnorm(mu1, mean(pooled), 1000*sd(pooled)) * dnorm(mu2, mean(pooled), 1000*sd(pooled)) * dexp(sig1, rate=0.1) * dexp(sig2, 0.1)
}
## Posterior
posterior <- function(parameters) {likelihood(parameters) * prior(parameters)}

## Sampling of posterior by MCMC
# parameters, starting values
pmu1 = 100; psig1 = 10; pmu2 = 100; psig2 = 10
parameters <- c(pmu1, psig1, pmu2, psig2)

#this is the MCMC /w Metropolis method
n.iter <- 10000
results <- matrix(0, nrow=n.iter, ncol=4)
results[1, ] <- parameters
for (iteration in 2:n.iter){
  candidate <- parameters + rnorm(4, sd=0.5)
  ratio <- posterior(candidate)/posterior(parameters)
  if (runif(1) < ratio) parameters <- candidate # Metropolis modification
  results[iteration, ] <- parameters
}

## Don't forget to burn-in the first x values to minimize the impact of the starting values
# burn-in
results_all <- results
results <- results[500:n.iter,]

## Results
mu1 <- results[,1]
sig1 <- results[,2]
mu2 <- results[,3]
sig2 <- results[,4]

hist(mu1 - mu2)
hist(sig1-sig2)

hdiff1 <- ggplot()+
  geom_histogram(aes(mu1-mu2),fill='blue')+
  theme_bw()
hdiff2 <- ggplot()+
  geom_density(aes(mu1-mu2),color='blue',fill='gold')+theme_bw()
#hdiff1
#hdiff2
multiplot(hs1,hs2,hdiff1,hdiff2,cols = 2)

## Answer questions with uncertainty
mean(mu1 - mu2 < 0)
mean(mu2 - mu1 > 1.0)
## Cumulative Density Function
qplot(sort(mu1-mu2),cumsum((sort(mu1-mu2))<0)/length(mu1))
## CDF without normalizing
qplot(sort(mu1-mu2),cumsum((sort(mu1-mu2))<0))

## Misc
# Plot delta (sorted)
qplot(1:length(mu1),sort(mu1 - mu2))
# Plot delta (unsorted)
qplot(1:length(mu1),(mu1 - mu2))


###################################################################
### RJAGS version ###
#####################

### ATTENTION: Requires installation of JAGS - Just Another Gibbs Sampler - first, can be downloaded from sourceforge 
### https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/
### If not installed, installing BEST and rjags still works, just loading the rjags package via require(rjags)/library(rjags) fails
library(BEST)
library(rjags)
#remove.packages("BEST")
#remove.packages("rjags")

## JAGS model definition - looks strange but works
model.str <- 'model {
    for (i in 1:Ntotal) {
y[i] ~ dt(mu[x[i]], tau[x[i]], nu)
}
for (j in 1:2) {
mu[j] ~ dnorm(mu_pooled, tau_pooled)
tau[j] <- 1 / pow(sigma[j], 2)
sigma[j] ~ dunif(sigma_low, sigma_high)
}
nu <- nu_minus_one + 1
nu_minus_one ~ dexp(1 / 29)
}'

# Indicator variable
x <- c(rep(1, length(sample.1)), rep(2, length(sample.2)))

cpd.model <- jags.model(textConnection(model.str),
                        data=list(y=pooled,
                                  x=x,
                                  mu_pooled=mean(pooled),
                                  tau_pooled=1/(1000 * sd(pooled))^2,
                                  sigma_low=sd(pooled) / 1000,
                                  sigma_high=sd(pooled) * 1000,
                                  Ntotal=length(pooled)))
update(cpd.model, 1000)
chain <- coda.samples(model = cpd.model, n.iter = 10000,
                      variable.names = c('mu', 'sigma'))
rchain <- as.matrix(chain)
#hist(rchain[, 'mu[1]'] - rchain[, 'mu[2]'])

## Probabilities for mu1 > mu2 and mu1+5 > mu2: 
mean(rchain[, 'mu[1]'] - rchain[, 'mu[2]'] < 0)
mean(rchain[, 'mu[2]'] - rchain[, 'mu[1]'] > 5)

## Plots for mean values
hdiff3 <- ggplot()+
  geom_histogram(aes(rchain[, 'mu[1]'] - rchain[, 'mu[2]']),fill='blue')+
  theme_bw()
hdiff4 <- ggplot()+
  geom_density(aes(rchain[, 'mu[1]'] - rchain[, 'mu[2]']),color='blue',fill='gold')+theme_bw()
#hdiff4
multiplot(hs1,hs2,hdiff3,hdiff4,cols = 2)
multiplot(hdiff1,hdiff2,hdiff3,hdiff4,cols = 2)

## Comparison to t-test
# t-test significant
t.test(sample.1,sample.2)
t.test(sample.1,sample.2-1)
mean(rchain[, 'mu[2]'] - rchain[, 'mu[1]'] > 1)
# t-test not significant
t.test(sample.1,sample.2-2)
mean(rchain[, 'mu[2]'] - rchain[, 'mu[1]'] > 2)

### Statistical Power etc.
library(pwr)
pwr.t.test(n=100,d=0.2,sig.level = 0.05)
pwr.t.test(n=1000,d=0.2,sig.level = 0.05)
pwr.t.test(n=200,d=0.2,sig.level = 0.05)
pwr.t.test(n=300,d=0.2,sig.level = 0.05)
pwr.t.test(n=400,d=0.2,sig.level = 0.05)
pwr.t2n.test(n1=300,n2=500,d=0.2,sig.level = 0.05)


### Investigate variance / standard deviations
hmu12 <- ggplot()+
  geom_histogram(aes(mu1),fill='blue',bins = 100)+
  geom_histogram(aes(mu2),fill='red',bins = 100)+
  xlab('Estimated means')+
  theme_bw()
hsig12 <- ggplot()+
  geom_histogram(aes(sig1),fill='blue',bins = 100)+
  geom_histogram(aes(sig2),fill='red',bins = 100)+
  xlab('Estimated standard deviations')+
  theme_bw()

multiplot(hmu12,hsig12,cols = 1)


## Illustrate burn-in effect

## Bootstrapping with the boot package
library(boot)
## TODO
