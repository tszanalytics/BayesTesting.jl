### functions for posterior odds and Bayes p-value
#
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
    pdr_pval(theta_draws)

Input:

        theta_draws = MC sample from posterior for theta

Optional argument:

        h0 = value in hull hypothesis (default = 0).

Returns:

        posterior density ratio, one-tailed and two-tailed 'p-value' of evidence against h0.
        one-tailed 'p-value' = minimum[Prob(theta < 0), Prob(theta>0)]
"""
# x = MC sample, h0 = null hypothesis value
# usage: odds, p_val, pval_2 = post_odds_pval(x, h0 = -2.3)
function pdr_pval(x; h0 = 0.0)
    d = kde(x)
    p,ind = findmax(d.density)
    ph0 = ifelse(pdf(d,h0) == 0.0,0.0000001,pdf(d,h0))
    odds = pdf(d,d.x[ind])/ph0
    p_val = length(x[x .<= h0])
    p_val2 = length(x[x .>= h0])
    if p_val <= p_val2
        p_value = p_val/length(x)
    else
        p_value = p_val2/length(x)
    end
    p_value_2tail = 2*p_value
    return odds, p_value, p_value_2tail
end



"""
    bayespval(theta_draws)

Input:

        theta_draws = MC sample from posterior for theta

Optional argument:

        h0 = value in hull hypothesis (default = 0).

Returns:

        posterior 'p-value' of evidence against h0.
        pval = one tailed value minimum[(Prob(theta <0), Prob(theta>0)]
        pval2 = two-tailed value (2*pval)
"""
function bayespval(mcs; h0=0.0)
  p_val = length(mcs[mcs .<= h0])/length(mcs)
  if p_val == 0.0
    p_val = 0.00001
  end
  if p_val > 0.5
    p_val = (1.0 - p_val)
  end
  p_val2 = 2.0*p_val
  return p_val, p_val2
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
