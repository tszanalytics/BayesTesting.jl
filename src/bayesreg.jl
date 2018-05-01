
function bayesreg(y,x)
  n = length(y)
  k = size(x,2)
  X = [ones(n) x]
  # uninformative analytical
  iXX = inv(X'*X)
  bhat = iXX*X'*y  # bhat = X \ y
  s2hat = (y - X*bhat)'*(y - X*bhat)/n
  Vb = s2hat.*iXX
  seb = sqrt.(diag(Vb))

  println("linear model results:")
  println("coeff estimates ", bhat)
  println("standard errors ",seb)
  println("s^2 (eqn. variance) = ",s2hat)

  # compute R^2
  tss = sum((y .- mean(y)).^2)
  R2 = 1 - s2hat*n/tss
  println("Rsquared = ",R2)
  yhat = X*bhat;

  return bhat, seb, s2hat, Vb, R2, yhat
end


function bayespval2(mcs; h0=0.0)
  sort!(mcs)
  p = 0
#  pvals = 0.0
  for i = 1:length(mcs)
    if mcs[i] <= h0
#      pvals += mcs[i]
      p += 1
    end
  end

  if p == 0
    pval = 0.0001
  elseif p == length(mcs)
    pval = 0.0001
  else
    pval = 2*p/length(mcs)
  end

  if pval >= 1.0
    pval = 2*(1.0 - pval/2)
  end
  return pval
end


# conjugate NIG prior regression:
"""
  bayesregNIG(y,x)

conjugate NIG prior linear regression function: Compute coeff. estimates, s.e's, equation σ, Rsquared

usage:

    bhat, seb, s2hat, Vb, R2, yhat = bayesregNIG(y,x; b0=0.,B0=100000.,a0=0.00001,d0=0.00001)
"""
function bayesregNIG(y,x; b0=0.,B0=100000.,a0=0.00001,d0=0.00001)
  n = length(y)
  X = [ones(n) x]
  k = size(X,2)

  # informative analytical- Greenberg, ch. 4, (4.5), (4.6)
  iB0 = inv(B0) # B0\1  # better way to invert computationally
  B1 = inv(X'*X + iB0)
  b1 = B1*(X'*y + iB0*b0)  # bhat = X \ y
  a1 = a0 + n - k
  postssr = (y - X*b1)'*(y - X*b1)
  d1 = d0 + postssr
  s2hat = d1/a1
  Vb = s2hat.*B1
  seb = sqrt.(diag(Vb))

  println("linear model posterior results:")
  println("coeff estimates ", b1)
  println("standard errors ",seb)
  println("s^2 (eqn. variance) = ",s2hat)

  # compute R^2
  tss = sum((y .- mean(y)).^2)
  R2 = 1 - s2hat*n/tss
  println("Rsquared = ",R2)
  yhat = X*b1

  return b1, seb, s2hat, Vb, R2, yhat
end

"""
  blinreg(y,x)

linear regression function: Compute coeff. estimates, s.e's, equation σ, Rsquared

usage:

  b, seb, s, R2 = linreg(y,x)
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
  R2 = 1 - s2*n/tss
  println("Rsquared = ",round(R2,3))
  s = sqrt(s2)
  return b, seb, s, R2;
end
