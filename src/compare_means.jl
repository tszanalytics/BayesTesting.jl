"""
Comparison of means: plot posterior densities and compute odds against difference of 0.

Returns:
        plt,diff_mean,qs,tst

        plt = plot of densities
        diff_mean = MC sample from posterior of difference in means
        qs = HPD quantile intervals from diff_mean
        tst = posterior odds against zero difference and area in tail to one side of zero.

Enter either:
        compare_mean(x,y)
        x, y are two samples to compare_means

        compare_means(m1, m2, s1, s2, n1, n2)
        m1, m2, are sample means, s1, s2, are sample SDs, n1, n2 are sample number of obs.

Optional arguments:
        M = MC sample size (larger for improvement numerical accuracy)
        lbl = ["mu 1 label" "mu 1 label"],  use lbl = ["" ""] to eliminate labels
        lgd = :topright (choose legend position, :topleft to move to top left of figure, etc.)
"""
function compare_means(m1, m2, s1, s2, n1, n2; M = 1000000,  lbl=["mu 1" "mu 2"],lgd = :topright)
    draws_m1 = m1 .+ (s1/sqrt(n1)).*rand(TDist(n1-1),M)
    draws_m2 = m2 .+ (s2/sqrt(n2)).*rand(TDist(n2-1),M)
    diff_mean = draws_m1 - draws_m2
# plot(draws_m1,st=:histogram,normalize=true,bins=100)
    l = @layout([a b])
    plt1 = plot(draws_m1,st=:density,fill=(0,0.4,:red),alpha=0.4,label=lbl[1],legend=lgd,title="Posteriors for means")
    plot!(draws_m2,st=:density,fill=(0,0.4,:blue),alpha=0.4,label=lbl[2])
    plt2 = plot(diff_mean,st=:density,fill=(0,0.4,:green),alpha=0.4,label="",title="Posterior difference")
    vline!([0.0],color=:black,label="")
    plt3 = plot(plt1, plt2, layout=l)
    qs = quantile(diff_mean,[0.005 0.025 0.05 0.5 0.95 0.975 0.995])
    @show(qs)
    tst = [mcodds(diff_mean) bayespval(diff_mean)]
    @show(tst)
    return plt3,diff_mean,qs,tst
end

# function overload for using two samples directly
function compare_means(x,y; M = 1000000, lbl=["mu 1" "mu 2"],lgd = :topright)
    m1 = mean(x); m2 = mean(y)
    s1 = std(x); s2 = std(y)
    n1 = length(x); n2 = length(y)
    draws_m1 = m1 .+ (s1/sqrt(n1)).*rand(TDist(n1-1),M)
    draws_m2 = m2 .+ (s2/sqrt(n2)).*rand(TDist(n2-1),M)
    diff_mean = draws_m1 - draws_m2
# plot(draws_m1,st=:histogram,normalize=true,bins=100)
    l = @layout([a b])
    plt1 = plot(draws_m1,st=:density,fill=(0,0.4,:red),alpha=0.4,label=lbl[1],legend=lgd,title="Posteriors for means")
    plot!(draws_m2,st=:density,fill=(0,0.4,:blue),alpha=0.4,label=lbl[2])
    plt2 = plot(diff_mean,st=:density,fill=(0,0.4,:green),alpha=0.4,label="",title="Posterior difference")
    vline!([0.0],color=:black,label="")
    plt3 = plot(plt1, plt2, layout=l)
    qs = quantile(diff_mean,[0.005 0.025 0.05 0.5 0.95 0.975 0.995])
    @show(qs)
    tst = [mcodds(diff_mean) bayespval(diff_mean)]
    @show(tst)
    return plt3,diff_mean,qs,tst
end

"""
 2x2 contigency table (categorical data) comparison of proportions
 """
function compare_proportions(s1, s2, n1, n2; a =1, b = 1,M = 1000000, lbl=["mu 1" "mu 2"],lgd = :topright)
    draws_m1 = rand(Beta((s1+a),(n1-s1+b)),M)
    draws_m2 = rand(Beta((s2+a),(n2-s2+b)),M)
    diff_mean = draws_m1 - draws_m2
# plot(draws_m1,st=:histogram,normalize=true,bins=100)
    l = @layout([a b])
    plt1 = plot(draws_m1,st=:density,fill=(0,0.4,:red),alpha=0.4,label=lbl[1],legend=lgd,title="Posteriors for proportions",ylab="Density",xlab="Proportion")
    plot!(draws_m2,st=:density,fill=(0,0.4,:blue),alpha=0.4,label=lbl[2])
    plt2 = plot(diff_mean,st=:density,fill=(0,0.4,:green),alpha=0.4,label="",title="Posterior difference",xlab="Difference in proportion")
    vline!([0.0],color=:black,label="",linewidth=2)
    plt3 = plot(plt1, plt2, layout=l)
    qs = quantile(diff_mean,[0.005 0.025 0.05 0.5 0.95 0.975 0.995])
    @show(qs)
    tst = [mcodds(diff_mean) bayespval(diff_mean)]
    @show(tst)
    return plt3,diff_mean,qs,tst
end

function compare_proportions(x, y; a =1, b = 1,M = 1000000, lbl=["mu 1" "mu 2"],lgd = :topright)
    n1 = length(x); s1 = sum(x)
    n2 = length(y); s2 = sum(y)
    draws_m1 = rand(Beta((s1+a),(n1-s1+b)),M)
    draws_m2 = rand(Beta((s2+a),(n2-s2+b)),M)
    diff_mean = draws_m1 - draws_m2
# plot(draws_m1,st=:histogram,normalize=true,bins=100)
    l = @layout([a b])
    plt1 = plot(draws_m1,st=:density,fill=(0,0.4,:blue),alpha=0.4,label=lbl[1],legend=lgd,title="Posteriors for proportions")
    plot!(draws_m2,st=:density,fill=(0,0.4,:red),alpha=0.4,label=lbl[2])
    plt2 = plot(diff_mean,st=:density,fill=(0,0.4,:green),alpha=0.4,label="",title="Posterior difference")
    vline!([0.0],color=:black,label="")
    plt3 = plot(plt1, plt2, layout=l)
    qs = quantile(diff_mean,[0.005 0.025 0.05 0.5 0.95 0.975 0.995])
    @show(qs)
    tst = [mcodds(diff_mean) bayespval(diff_mean)]
    @show(tst)
    return plt3,diff_mean,qs,tst
end

"""
Function to draw from Beta(s+a, n-s+b) posterior

     Default prior is uniform: a = b = 1
"""
function beta_posterior(s,n; a=1, b=1, M=100000)
    theta_draws = rand(Beta((s+a),(n-s+b)),M)
    return theta_draws
end

"""
Function to update Beta(s+a, n-s+b) posterior given two samples

     Default prior is uniform: a = b = 1 for each sample distribution
"""
function beta_update(s1,n1,s2,n2; a1=1,b1=1,a2=1,b2=1,M=100000)
    theta_draws = rand(Beta((s1+s2+a1+a2),(n1+n2-s1-s2+b1+b2)),M)
    return theta_draws
end


"""
equiv_test = Equivalence test from MC sample of difference in ATE

    mc_sample = posterior simulation sample of difference in treatments

    interval = 1/2 tolerance interval  = distance from 0.0 to one side of interval

    see TADS_equiv.jl for example of use.
"""
function equiv_test(mc_sample,interval;h0=0.0)
  p1 = p2 = p3 = p4 = p5 = 0.0
  for th1 in mc_sample
      if th1 < (h0-3*interval)
          p1 += 1
      elseif th1 < (h0-interval)
          p2 += 1
      elseif th1 < (h0+interval)
          p3 += 1
      elseif th1 < (h0+3*interval)
          p4 += 1
      else
          p5 += 1
      end
  end
  ps = [p1; p2; p3; p4; p5]
  normps = ps./sum(p1 + p2 + p3 + p4 + p5)
  x = [-4*interval; -2*interval; 0.0; 2*interval; 4*interval]
  plt = plot(mc_sample,st=:density,label="",ylab="",alpha=0.4,fill=(0,0.4,:green),ylim=(0.0,0.4))
  plot!(twinx(),x,normps, st=:bar, color=:grey,alpha=0.4,bar_width=2*interval,label="",ylim=(0.0,1.0),ylab="Probability",title = "Equivalence interval probabilities") #,axis=:right
  oddsf35 = normps[3]/normps[5]
  oddsf34 = normps[3]/normps[4]
  oddsf32 = normps[3]/normps[2]
  oddsf31 = normps[3]/normps[1]
  eqodds = [oddsf31 oddsf32 oddsf34  oddsf35]
  @show(eqodds)
  @show(normps)
  return eqodds, normps, plt
end
