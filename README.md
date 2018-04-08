# BayesTesting.jl
**Bayesian Hypothesis Testing without Tears**

Objective Bayesian hypothesis testing that does not suffer from the problems inherent in the standard approach and so work well in practice:
1.	The JLB paradox does not arise.
2.	Any prior can be employed, including uninformative priors, so the same prior employed for inference can be used for testing.
3.	In standard problems when the posterior distribution matches (numerically) the frequentist sampling distribution or likelihood, there is a one-to-one correspondence with the frequentist test.
4.	Provides posterior odds against the null hypothesis that are easy to interpret (unlike *p*-values), do not violate the likelihood principle, and result from minimizing a linear combination of type I and II errors rather than fixing the type I error before testing (as in Neyman-Pearson significance testing).







[![Build Status](https://travis-ci.org/Jeff Mills/BayesTesting.jl.svg?branch=master)](https://travis-ci.org/Jeff Mills/BayesTesting.jl)

[![Coverage Status](https://coveralls.io/repos/Jeff Mills/BayesTesting.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Jeff Mills/BayesTesting.jl?branch=master)

[![codecov.io](http://codecov.io/github/Jeff Mills/BayesTesting.jl/coverage.svg?branch=master)](http://codecov.io/github/Jeff Mills/BayesTesting.jl?branch=master)
