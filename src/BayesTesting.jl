module BayesTesting

using Distributions, PlotRecipes, KernelDensity, LinearAlgebra, Statistics

export
    update_mean,
    marginal_posterior_mu,
    mcodds,
    todds,
    bayespval,

    blinreg,
    bayesregNIG,
    bayesreg,
    gsreg,
    compare_means,
    compare_proportions,

    beta_posterior,
    beta_update,
    equiv_test,

    bayes_ttest,
    bayes_correlation



include("bayes_meta_analysis.jl")
include("bayes_odds_pvalue.jl")
include("bayesreg.jl")
include("gsreg.jl")
include("compare_means.jl")
include("bayes_ttest.jl")
include("bayes_correlation.jl")

end
