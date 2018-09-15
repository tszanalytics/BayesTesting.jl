using Distributions, KernelDensity, PlotRecipes, Statistics






"""
    bayes_ttest(x::Vector{Real}; h0::Real=0.0, iter::Real=10e4)

Perform a Bayesian t-test.

# Examples
```julia
x = randn(50) .+ 0.2
todds, pval, qs = bayes_ttest(x,h0=0.0)
@show(todds, pval)
@show(qs)
```
"""
function bayes_ttest(x::Vector{Real}; h0::Real=0.0, iter::Real=10e4)
	n = length(x)
	v = n - 1
	mu = mean(x)
	se = std(x)/sqrt(n)
	vm = var(x)/n
	tsq = ((mu - h0)^2)/vm
	todds = (1.0 + tsq/v)^(0.5*(v+1))
	pval1 = cdf(TDist(v),sqrt(tsq))
	pval = 2.0*minimum([pval1,(1.0 - pval1)])
	ts = se.*rand(TDist(v), Int(iter)) .+ mu
	ps = [0.005; 0.025; 0.5; 0.975; 0.95]
	q = quantile(ts,[0.005,0.025,0.5,0.975,0.95])
	qs = [ps q]

	return todds, pval, qs
end





# Placeholder for the plot() method when types will be available
function plot_bayes_ttest(mu, ts, h0)
	plt = plot(ts,st=:density, fill=(0,0.4,:blue),alpha=0.4,label="posterior")
	vline!([mu h0],linecolor = [:green :black],label=["mean" "H0"])
    return plt
end
