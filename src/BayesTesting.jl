module BayesTesting

export update_mean, marginal_posterior_mu, mcodds,
        todds, bayespval, blinreg, bayesregNIG, bayesreg, gsreg
# package code goes here
include("bayes_meta_analysis.jl")
include("bayes_odds_pvalue.jl")
include("bayesreg.jl")
include("gsreg.jl")

end # module
