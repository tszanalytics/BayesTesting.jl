using Distributions, PlotRecipes





"""
    bayes_correlation(x, y; h0=0.0)

Perform a Bayesian correlation test.

# Examples
```julia
x = randn(50)
y = 0.5 .+ 0.3.*x .+ randn(50)

todds, pval, qs = bayes_correlation(x, y)
```
"""
function bayes_correlation(x::Vector{Real}, y::Vector{Real}; h::Real0=0.0, iter::Real=10e4)
    n = length(y)
    v = n - 2
    X = [ones(n) x]
    b = inv(X'X)*X'y
    residuals = y .- X*b

    s2 = sum(residuals.^2)/v
    vb = s2.*inv(transpose(X)*X)
    seb = sqrt(vb[2,2])
    that = abs(b[2]/seb)
    todds = (1.0 + (that^2)/v)^(0.5*(v+1))
    pval = 2.0*minimum([cdf(TDist(v),that) (1.0 - cdf(TDist(v),that))])
    ts = (seb*std(x)/std(y)).*rand(TDist(v), Int(iter)) .+ (b[2]*std(x)/std(y))
	ps = [0.005; 0.025; 0.5; 0.975; 0.95]
	q = quantile(ts,[0.005,0.025,0.5,0.975,0.95])
	qs = [ps q]

	return todds, pval, qs
end






# Placeholder for the plot() method when types will be available
function plot_bayes_correlation(x, y, ts, h0)
    plt = plot(ts,st=:density, label="posterior", fill=(0,0.4,:blue), alpha=0.4, title="Posterior for Correlation Coefficient")
	vline!([cor(y,x) h0], linecolor = [:green :black], label=["mean" "H0"])
    return plt
end
