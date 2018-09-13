using BayesTesting
using Base.Test

# write your own tests here
@testset "Basic tests" begin
    Random.seed(1235)
    @test round(sum(marginal_posterior_mu(0,1,10,M=1000)),3) == 42.416

    ms = [10 8 5]
    sds = [3 2 2]
    ns = [100 100 100]
    @test standardize_scores(ms,sds,ns) == ([1.0 0.8 0.5], [0.03 0.02 0.02])

    @test round(sum(update_mean(0.0,1.0,1.0,1.0,20,20)),3) == 41.731

end
