# BayesTesting.jl
**Bayesian Hypothesis Testing without Tears**

Objective Bayesian hypothesis testing that does not suffer from the problems inherent in the standard approach and so work well in practice:
1.	The JLB paradox does not arise.
2.	Any prior can be employed, including uninformative priors, so the same prior employed for inference can be used for testing.
3.	In standard problems when the posterior distribution matches (numerically) the frequentist sampling distribution or likelihood, there is a one-to-one correspondence with the frequentist test.
4.	Provides posterior odds against the null hypothesis that are easy to interpret (unlike *p*-values), do not violate the likelihood principle, and result from minimizing a linear combination of type I and II errors rather than fixing the type I error before testing (as in Neyman-Pearson significance testing).

**Functions currently available (package is under development)**

**Hypothesis testing:**

Optional parameter in following functions: h0= value in hull hypothesis (default is h0 = 0)

todds(theta_hat,theta_hat_se,v) = returns Student-t posterior odds for theta
todds(theta_draws,v) = returns Student-t posterior odds for theta given MC sample for theta

mcodds(theta_draws) = returns posterior odds given MC sample for theta (any distribution).

bayespval(theta_draws) = returns Bayesian p-value (tail area) give MC sample for theta

**Posterior inference:**

update_mean(m1,m0,s1,s0,n1,n0) = For Gaussian posterior sample 1 (or prior) with mean = m0, sd = s0, number of obs. =n0, and Gaussian likelihood or posterior for sample 2 with mean = m1, SD = s1, number of obs. = n1, returns tuple of combined sample posterior mean = m2, SD = s2, number of obs. = n2

marginal_posterior_mu(m,s, n, M) = return M draws from Student-t marginal posterior density with mean = m, SD = s, number of obs. = n.  M is an optional argument (default is M = 10000).

blinreg(y,X) = estimate a linear model y=Xβ+u (define X to contain vector of ones for an intercept)


Example 1: Testing if a sample mean equals zero

```
using BayesTesting

n = 50
x = randn(n)    # generate data

v = n-1           # degrees of freedom
t_hat = sqrt(n)*mean(x)/std(x)  # t-statistic
todds(t_hat,v)   # posterior odds vs null = 0

# Result: todds(t_hat, v) = 2.827
```

**Installation**
Currently unregistered, to install use Pkg.clone() with the repository url:
`Pkg.clone("git@github.com:tszanalytics/BayesTesting.jl.git")`


