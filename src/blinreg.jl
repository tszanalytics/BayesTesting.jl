
"""
linear regression function
  Compute coeff. estimates, s.e's, equation σ, Rsquared

usage:

    b, seb, s, R2 = linreg(y,x)
    y = dependent variable vector
    x = matrix of independent variables (no intercept in x)
"""
function blinreg(y,x)
  # add intercept
  n = length(y)
  X = [ones(n) x]
  b = (X'*X) \ X'*y
  resids = y - X*b
  RSS = sum(resids.^2)
  s2 = RSS/n

  covb = s2.*inv(X'X)
  seb = sqrt.(diag(covb))
  k = length(b)
  ### Correct formulas for AIC and BIC
  BIC = n*log(s2) + k*log(n)
  AIC = 2*k + n*log(s2)
  odds = zeros(k)
  pvals = zeros(k)
  that = abs.(b)./seb
  for i in 1:k
    odds[i] = (1.0 + (that[i]^2)/(n-k))^((n-k+1)/2.0)
    pvals[i] = 1.0 - cdf(TDist(n-k), that[i])
  end
  println("coeffs = ", round.(b, digits=3))
  println(" s.e's = ", round.(seb, digits=3))
  println(" odds = ", round.(odds, digits=3))
  println("p-values = ", round.(pvals, digits=4))
  println("s^2 (eqn. variance) = ",round(s2, digits=6))

  # compute R^2
  tss = sum((y .- mean(y)).^2)

  R2 = 1 - s2*n/tss
  println("Rsquared = ",round(R2, digits=3))
  println("AIC = ", round(AIC, digits = 2), "  BIC = ", round(BIC, digits = 2),)
  s = sqrt(s2)
  return b, seb, odds, pvals, s, R2, RSS, AIC, BIC
end

## example of use:
#x = randn(20)
#y = 1 + 1.*x .+ randn(20)
#b, seb, odds, pvals, s, R2, RSS, AIC, BIC = blinreg(y,x)
