# using Distributions, PlotRecipes

# Order of variables in function does not matter (since cor(y,x) = cor(x,y))
function correlation_ttest(y,x; h0=0.0)
    n = length(y)
    v = n - 2
    X = [ones(n) x]
    b = inv(X'X)*X'y
    res = y .- X*b
    s2 = sum(res.^2)/v
    vb = s2.*inv(transpose(X)*X)
    seb = sqrt(vb[2,2])
    that = abs(b[2]/seb)
    todds = (1.0 + (that^2)/v)^(0.5*(v+1))
    pval = 2.0*minimum([cdf(TDist(v),that) (1.0 - cdf(TDist(v),that))])
    ts = (seb*std(x)/std(y)).*rand(TDist(v),1000000) .+ (b[2]*std(x)/std(y))
	ps = [0.005; 0.025; 0.5; 0.975; 0.95]
	q = quantile(ts,[0.005,0.025,0.5,0.975,0.95])
	qs = [ps q]
	plt = plot(ts,st=:density, label="posterior", fill=(0,0.4,:blue),alpha=0.4, title="Posterior for Correlation Coefficient")
	vline!([cor(y,x) h0],linecolor = [:green :black],label=["mean" "H0"])
	return todds, pval, qs, plt
end


# Example of use (vary coefficient on x to increase/reduce correlation)
#x = randn(50)
#y = 0.5 .+ 0.3.*x .+ randn(50)
#cor(x,y)
#results = correlation_ttest(x,y)
#results[1]
#results[2]
#results[3]
#results[4]
# savefig("trash.png")
