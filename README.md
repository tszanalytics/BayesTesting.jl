# BayesTesting.jl
**Bayesian Hypothesis Testing without Tears**

Objective Bayesian hypothesis testing that does not suffer from the problems inherent in the standard approach and so work well in practice:
1.	The JLB paradox does not arise.
2.	Any prior can be employed, including uninformative priors, so the same prior employed for inference can be used for testing.
3.	In standard problems when the posterior distribution matches (numerically) the frequentist sampling distribution or likelihood, there is a one-to-one correspondence with the frequentist test.
4.	Provides posterior odds against the null hypothesis that are easy to interpret (unlike *p*-values), do not violate the likelihood principle, and result from minimizing a linear combination of type I and II errors rather than fixing the type I error before testing (as in Neyman-Pearson significance testing).

**Installation**
Currently unregistered, to install use Pkg.clone() with the repository url:

`Pkg.clone("https://github.com/tszanalytics/BayesTesting.jl")`


**Functions currently available (package is under development)**

**Hypothesis testing:**

Optional parameter in following functions: h0= value in hull hypothesis (default is h0 = 0)

**todds(theta_hat,theta_hat_se,v)** = returns Student-t posterior odds for theta

**mcodds(theta_draws)** = returns posterior odds given MC sample for theta (any distribution).

**bayespval(theta_draws)** = returns Bayesian p-value (tail area) give MC sample for theta

**Posterior inference:**

**update_mean(m1,m0,s1,s0,n1,n0)** = For Gaussian posterior sample 1 (or prior) with mean = m0, sd = s0, number of obs. =n0, and Gaussian likelihood or posterior for sample 2 with mean = m1, SD = s1, number of obs. = n1, returns tuple of combined sample posterior mean = m2, SD = s2, number of obs. = n2

**marginal_posterior_mu(m,s, n, M)** = return M draws from Student-t marginal posterior density with mean = m, SD = s, number of obs. = n.  M is an optional argument (default is M = 10000).

**blinreg(y,X)** = estimate a linear model y=Xβ+u (define X to contain vector of ones for an intercept)

**gsreg(y,X)** = Gibbs sampler for linear regression with default uninformative prior, X must contain 
	vector of ones to include intercept.
                Optional parameters:
                tau = precision starting value (default = 1.0)
                M = MCMC sample size (default = 10,000)

**gsreg(y,X, M=m, tau=t, b0=priorb, iB0 = invpriorcovb , d0=b, a0=a)** = Gibbs sampler with NIG prior.
          Note: iB0 = prior precision matrix = inv(prior variance matrix)
                b0 must be a column vector, 
                a0 and b0 are prior parameters for tau ~ Gamma(a,b)


**Example 1: Testing if a sample mean equals zero**

```
using BayesTesting
srand(1235)             # generate psuedo-data, n obs.
n = 50
x = randn(n)

v = n-1                 # degrees of freedom
mu_hat = mean(x)        # sample mean
se_mu = std(x)/sqrt(v)  # sample standard error of mean
todds(mu_hat,se_mu,v)   # posterior odds vs. zero

# Result: todds(mu_hat, se_mu, v, h0=0) = 1.016  => 1:1 odds against the null.


# with a nonzero mean - change the data generating process for x above to:
x = 0.5 + randn(n)
# Resulting posterior odds: todds(mu_hat, se_mu, v, h0=0) = 110.50  => 110:1 odds against the null
```

**More detailed help and examples in:** BayesTesting.jl_docs_2018.pdf


