module BayesTesting

export standardize_scores, update_mean, marginal_posterior_mu, mcodds, todds, bayespval
# package code goes here
include("bayes_meta_analysis.jl")
include("bayes_odds_pvalue.jl")

end # module
