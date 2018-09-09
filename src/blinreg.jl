
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

  s2 = sum(resids.^2)/n
  covb = s2.*inv(X'X)
  seb = sqrt.(diag(covb))
  k = length(b)
  odds = zeros(k)
  for i in 1:k
    odds[i] = todds(b[i],seb[i],(n-k))
  end
  println("coeffs = ",round.(b,3))
  println(" s.e's = ",round.(seb,3))
  println(" odds = ",round.(odds,3))
  println("s^2 (eqn. variance) = ",round(s2,6))

  # compute R^2
  tss = sum((y .- mean(y)).^2)
  R2 = 1.0 - s2*n/tss
  println("Rsquared = ",round(R2,3))
  s = sqrt(s2)
  return b, seb, odds, s, R2
end

## example of use:
#x = randn(20)
#y = 1 + 1.*x .+ randn(20)
#b, seb, odds, s,R2 = blinreg(y,x)
