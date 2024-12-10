
library(Bhat) # R library needed for calculation 
              # of profile LRCIs.

# We first record the values of the parameters 
# we need to calculate the profile LRCIs for the 
# binomial and beta parameters.  For the binomial 
# profile LRCI,  'k' must be greater than zero and 
# less than n.  For the beta profile LRCI, 'a' and
# 'b' must be greater than or equal to 1, 
n = 20 # sample size
k = 1 # number of successes
p.hat = k/n # MLE of binomial parameter 'p'
cl = 0.95   # Use cl = 0.956151315  to match the 
            # HPD intervals for the later example.
prior1 = 1 # uniform prior for beta 'a' parameter
prior2 = 1 # uniform prior for beta 'b' parameter
a = k + prior1  # beta parameter a
b = n - k +  prior2 # beta parameter b
mu = a/(a + b) # beta mean

# Now we calculate the MLE for the parameter 'p' in 
# the binomial distribution.  First, we create the 
# negative log-likelihood function to be optimized.
# Then we fit the distribution parameters and find
# the MLE.
binom.nll <- function(p) {
  return(-1*log(dbinom(k, size=n, prob=p)))
}

binom.fit <- nlminb(p.hat, binom.nll, lower=0, 
                    upper=1)
binom.MLE <- binom.fit$par

# Calculate LB and UB for 95% LRCI for the binomial 
# parameter 'p'.
plist=list(label="p",est=binom.MLE, low=0,upp=1)
binom.LRCI=plkhci(plist,binom.nll,"p",prob=cl, 
                  eps=1e-16)

# Show the calculated values of the MLE and profile 
# LRCI for the binomial parameter 'p'.
binom.MLE
round(binom.LRCI,6) #Print to 6 decimal places.

# Calculate the MLE for the parameter 'x' in the 
# beta distribution.  First, create the negative 
# log-likelihood function to be optimized.
beta.nll <- function(x) {
  return(-1*log(dbeta(x, shape1 = a, shape2 = b)))
}

beta.fit <- nlminb(p.hat, beta.nll, lower=0, upper=1)
beta.MLE <- beta.fit$par

# Calculate LB and UB for 95% LRCI for the beta 
# parameter 'x'.
plist=list(label="x",est=beta.MLE, low=0,upp=1)
beta.LRCI=plkhci(plist,beta.nll,"x",prob=cl, eps=1e-16)

# Show the calculated values of the MLE and profile 
# LRCI  for the beta parameter 'x'.
beta.MLE
round(beta.LRCI, 6) #Print to 6 decimal places.

# Now we calculate the HPD intervals for the 
# beta distribution.
library(HDInterval) # R library needed for calculation 
                    # of HPD intervals.
cl=0.95             # Use cl = 0.9432647 to match the 
                    # profile LRCIs from earlier.

# Calculate the HPD interval for the beta distribution.
round(hdi(qbeta, cl, shape1=a, shape2=b), 6)

