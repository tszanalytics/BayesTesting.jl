### functions for posterior odds and Bayes p-value
#
using KernelDensity

"""
    todds(theta_hat,theta_hat_se,v)

Input:

        v = degrees of freedom

        theta_hat, theta_hat_se
        = estimate and standard error for theta

Optional argument:

        h0 = value in hull hypothesis (default = 0).


Returns:

        t-posterior odds ratio of evidence against h0 (default = 0).
"""
function todds(theta_hat, theta_se, v; h0 = 0.0)
  # compute posterior odds assuming a Student-t distribution
  # mcs = MC sample
  # v = degrees of freedom, n - k
  # h0 = value under null hypothesis
  that = abs.(theta_hat - h0)/theta_se
  odds = (1.0 + (that^2)/v)^((v+1)/2.0)
end

"""
    mcodds(theta_draws)

Input:

        theta_draws = MC sample from posterior for theta

Optional argument:

        h0 = value in hull hypothesis (default = 0).

Returns:

        posterior odds ratio of evidence against h0.
"""
function mcodds(mcs; h0=0.0)
  # dependency: KernelDensity
  # compute posterior odds from nonparametric density
  # mcs = MC sample

### NEED catches for when odds too small/large - out of bounds of sample
  postden = KernelDensity.kde(mcs)
  ax = abs.(postden.x .- h0)
  indx = findall(ax .== minimum(ax))  # find all x that = 4 in x
  numodds = maximum(postden.density)
  denodds = postden.density[indx[1]]
  odds = numodds/denodds[1]
  if odds > 5000000.0
    odds = 5000000.0
  end
  return odds
end

# Example: if bs[:,2] is vector of MCMC draws for beta
# to test if beta = 1.00:
#mcodds(bs[:,2],h0=1.00)

"""
    bayespval(theta_draws)

Input:

        theta_draws = MC sample from posterior for theta

Optional argument:

        h0 = value in hull hypothesis (default = 0).

Returns:

        posterior 'p-value' of evidence against h0.
"""
function bayespval(mcs; h0=0.0)
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
    pval = 2.0*p/length(mcs)
  end

  if pval >= 1.0
    pval = 2.0*(1.0 - pval/2)
  end
  return pval
end

#= old deprecated todds assuming t
function todds(mcs,n,k; h0 = 0.0)
  # compute posterior odds assuming a Student-t distribution
  # mcs = MC sample
  # v = degrees of freedom, n - k
  # h0 = value under null hypothesis
  v = n-k
  unbiasedseb = std(mcs)*sqrt(v)/sqrt(n)
  that = abs.(mean(mcs) - h0)/unbiasedseb
  odds = (1 + (that^2)/v)^((v+1)/2)
end
=#
