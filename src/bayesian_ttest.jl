# using Distributions, KernelDensity, PlotRecipes, Statistics
function Bayesian_ttest(x; h0=0.0)
	n = length(x)
	v = n - 1
	mu = mean(x)
	se = std(x)/sqrt(n)
	vm = var(x)/n
	tsq = ((mu - h0)^2)/vm
	todds = (1.0 + tsq/v)^(0.5*(v+1))
	pval1 = cdf(TDist(v),sqrt(tsq))
	pval = 2.0*minimum([pval1,(1.0 - pval1)])
	ts = se.*rand(TDist(v),1000000) .+ mu
	ps = [0.005; 0.025; 0.5; 0.975; 0.95]
	q = quantile(ts,[0.005,0.025,0.5,0.975,0.95])
	qs = [ps q]
	plt = plot(ts,st=:density, fill=(0,0.4,:blue),alpha=0.4,label="posterior")
	vline!([mu h0],linecolor = [:green :black],label=["mean" "H0"])
	return todds, pval, qs, plt
end



# Example of use
#using Random
#Random.seed!(1235)
#x = randn(50) .+ 0.2
#t_odds, p_val, qs, plt = Bayesian_ttest(x,h0=0.0)
#@show(todds, pval)
#@show(qs)
#plt
# savefig("trash.png")
