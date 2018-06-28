# bayes_meta_analysis.jl
# functions for performing a Bayesian meta analysis
# Pkg.checkout("Plots","dev")
using Distributions  #,  Plots, StatPlots
# Normalize values for each experiment to allow comparison

function standardize_scores(means,stds,nobs)
  # Normalizing scores to [0,1] interval
  # Dividing by week 0 score

  # if no data, return all zeros
  for i in 1:length(nobs)
    if nobs[i] == 0
      nobs[i] = 1
    end
  end

  stdz_mn = means./means[1]
  # scale SD
  stdz_sd = stds./means[1]
  # standard error of scaled mean
  stdz_se = stdz_sd./sqrt.(nobs)
  return stdz_mn, stdz_se
end

### Cannot divide by std for inference!
#function standardize_scores_sddiv(means,stds,nobs)
  # Normalizing scores to [0,1] interval
  # Dividing by week 0 score
#  stdz_mn = (means[1] .- means) ./ stds
  # scale SD
  #stdz_sd = stds/means[1]
  # standard error of scaled mean
#  stdz_se = stdz_sd./sqrt(nobs)
#  return stdz_mn, stdz_se
#end



"""
        update_mean(m1,m0,s21,s22,n1,n0)

Compute updated (combined) posterior summary statistics given
two sets of posterior statistics assuming both Gaussian.

 - For uninformative prior, set m0 = s20 = n0 = 0.0
 - see Gelman et al. (2014), eqn. (3.7), p. 68

Input:

        m0 = posterior mean for first sample (or prior mean)

        m1 = posterior mean for second sample (or likelihood mean)

        s0 = posterior standard error for sample 1 mean (NOT SD of sample)

        s1 = posterior standard error for sample 1 mean (NOT SD of sample)

        n0 = number of obs. in sample 1

        n1 = number of obs. in sample 2

N.B.: s = std(x)/sqrt(n)

Returns:

        m2 = posterior mean of combined sample

        s2 = posterior standard error for combined sample

        n2 = n0 + n1
"""
function update_mean(m1,m0,s1,s0,n1,n0)
  s21 = s1^2
  s20 = s0^2
  n2 = n1 + n0
  m2 = (n1/n2)*m1 + (n0/n2)*m0
  v2 = n2 - 1.
  vs22 = (n1-1.)*s21 + (n0-1.)*s20 + ((n1*n0)/n2)*(m1 - m0)^2.
  s2 = sqrt(vs22/v2)
  return m2, s22, n2
end

"""
        marginal_posterior_mu(m,s,n, M)

  Compute M marginal Student-t posterior draws.

Input:

          posterior mean, m, standard error, s, and number of obs. n

          M = number of draws (optional: default M = 10000)

Returns:

          M draws for Student-t posterior density.
"""
# using Distributions
function marginal_posterior_mu(m,s,n; M = 10000)
  v = n - 1
  ts = m .+ s.*rand(TDist(v),M)
  return ts
end


"""
draw from posteriors for treatment and nontreatment
mn = mean, s2 = variance, n = no. of obs.
"""
function post_treat_notreat(mntreat,s2treat,ntreat, mnnotreat, s2notreat,nnotreat; M = 10000)
  # compute posteriors for treatment and nontreatment
  # group scores for each week/period
  wks = length(mntreat)
  drawst = zeros(M,wks)
  drawsnt = zeros(M,wks)
  for i in 1:wks
    drawst[:,i] = marg_post_mu(mntreat[i],s2treat[i],ntreat[i], M = M)
    drawsnt[:,i] = marg_post_mu(mnnotreat[i],s2notreat[i],nnotreat[i], M = M)
  end
  return drawst, drawsnt
end

"""
Compute posterior difference in treatment vs. nontreatment means
"""
function treat_pcbo_difference(drawst,drawsnt)
  difbywk = zeros(length(drawst[:,1]),length(drawst[1,:]))
  for i in 1:length(drawst[1,:])
    difbywk[:,i] = drawst[:,i] - drawsnt[:,i]
  end
  return difbywk
end

"""
Trajectory plot with 0.95 CIs given posterior means
"""
# change x axis labels to weeks instead of period
function traj_plot(means;overlay = false,nwks=3,ci=true)
  if (nwks == 3)
    z = [1 2 3]'   # for up to week 8 = 3 data points
  else
    z = [1 2 3 4]'   # for up to week 12 = 4 data points
  end

  mn = means
  lz = log.(z)
  Xloglin = hcat(ones(nwks),lz)
  b = Xloglin \ mn[1:nwks]
  s2 = sum((mn[1:nwks] - Xloglin*b).^2)/(nwks-1.)
  println("coeffs = ",b)
  println("s^2 (eqn. variance) = ",s2)

  # compute R^2
  tss = sum((mn[1:nwks]' - mean(mn[1:nwks])).^2)
  R2 = 1 - s2*2/tss
  println("Rsquared = ",R2)
  zs = 1.0:0.1:4.2
  score = b[1] + b[2].*log.(zs')

  #plot(1.0:0.1:4.2,score',size=(400,200),label="")
  mns = mn[1:nwks]
  #scatter!(z',mns',label="")

  s = sqrt(s2)
  tcrit5 = 4.3  # 5% crit. value from Student-t with 2 d.f.
  up95 = score .+ s*tcrit5
  lo95 = score .- s*tcrit5
  #  zz = 0.0:0.37:12.2
  zz = 1.0:0.1:4.2
  if (overlay == false)
    labels = [0,4,8,12]
#    zz = 0.0:0.37:12.2
 zz = 1.0:0.1:4.2
    plot(zz,score',size=(400,200),label="")
  else
    plot!(zz,score',size=(400,200),label="")
  end
  wks = [0 4 8 12]
  #scatter!(wks[1:nwks]',mns',label="")
  scatter!(z',mns',label="")

  if ci == true
    plot!(zz,up95',size=(400,200),label="",color=:green)
    plot!(zz,lo95',size=(400,200),label="",color=:green)
  end
#  hline!([1.0], color=:red,label="")
end


function pvals_by_week(m)
# p-value computation:
  df = csv_plot(m)
  x=df[:,1]   # scale horizontal axis to weeks
  wks = convert_to_12weeks(x)
  b,s,r2 = traj_reg_stats(m)
  score = df[:,2]   # predictive means
  bpval = zeros(length(score))
  for i in 1:length(score)
    #tvals = marg_post_mu(score[i],s^2,nwks)
    #bpval[i] = bayespval(tvals)
    t = score[i]/s
    bpval[i] = cdf(TDist(3),t)

  end
  zz = zeros(7,3)
  j = 1
  for i = 1:length(score)
    if (wks[i] == 0.0) || (wks[i] == 2.0) || (wks[i] == 4.0) || (wks[i] == 6.0) || (wks[i] == 8.0) || (wks[i] == 10.0) || (wks[i] == 12.0)
      println([ round(score[i],3),wks[i],bpval[i] ])
      zz[j,1], zz[j,2], zz[j,3] =  round(score[i],3), wks[i], bpval[i]
      j += 1
    end
  end
  return zz
end

function traj_cis(means;nwks=3)
  if (nwks == 3)
    z = [1 2 3]'   # for up to week 8 = 3 data points
  else
    z = [1 2 3 4]'   # for up to week 12 = 4 data points
  end

  mn = means
  lz = log.(z)
  Xloglin = hcat(ones(nwks),lz)
  b = Xloglin \ mn[1:nwks]
  s2 = sum((mn[1:nwks] - Xloglin*b).^2)/(nwks-1.)
  zs =  [1 2 3 4]'
  score = b[1] + b[2].*log.(zs')

  s = sqrt(s2)
  tcrit5 = 4.3  # 5% crit. value from Student-t with 2 d.f.
  up95 = score .+ s*tcrit5
  lo95 = score .- s*tcrit5
  return lo95, up95
end


function cond_post_mu(m,s2; M = 10000)
  # need to either evaluate pdf of t, or
  # draw from t and plot pseudo-sample
  ts = rand(Normal(m,sqrt(s2)),M)
  return ts
end


function traj_reg_stats(means; nwks=3)
  if (nwks == 3)
    z = [1 2 3]'   # for up to week 8 = 3 data points
  else
    z = [1 2 3 4]'   # for up to week 12 = 4 data points
  end

  mn = means
  lz = log.(z)
  Xloglin = hcat(ones(nwks),lz)
  b = Xloglin \ mn[1:nwks]
  s2 = sum((mn[1:nwks] - Xloglin*b).^2)/(nwks-1.)
  println("coeffs = ",b)
  println("s^2 (eqn. variance) = ",s2)

  # compute R^2
  tss = sum((mn[1:nwks]' - mean(mn[1:nwks])).^2)
  R2 = 1 - s2*2/tss
  println("Rsquared = ",R2)
  s = sqrt(s2)
  return b, s, R2
end
