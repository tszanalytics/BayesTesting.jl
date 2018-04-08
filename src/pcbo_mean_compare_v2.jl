# data from Dosing Extracted Data 2-12-16
# Walkup pbo
week = [0 4 8 12]
# nwks = length(week)
### use only up to week 8
nwks = 3

wpbomn = [19.6 16.0 13.6 12.6]
wpbosd = [3.9 4.1 5.2 6.3]
wnpbo = [76 76 76 76]  # all weeks

# walkup sert {Note: 1 = SSRI, 0 = SNRI (at end of name)
ws1mn = [18.8 14.2 11.2 9.8]
ws1sd = [3.9 4.0 5.0 6.2]
wns1 = [133 133 133 133] # all weeks

# Birmaher pbo
bpbomn = [14.9 11.4 8.6 9.3]
bpbosd = [3.5 3.7 5.1 4.8]
bnpbo = [37 37 37 37]

# Birmaher Flx
bs1mn = [15.6 9.8 7.8 7.1]
bs1sd = [3.5 5.4 5.4 5.9]
bns1 = [37 37 37 37]

# Wagner pbo
wgpbomn = [77.55 65.4 57.58 54.0]
wgpbosd = [19.34 21.44 12.15 25.8]
wgnpbo = [156 150 150 150]

# Wagner Parox
wgs1mn = [77.54 53.46 38.87 32.0]
wgs1sd = [22.34 21.69 21.29 21.1]
wgns1 = [163 159 159 159]

# Rynn pbo - no week 12 here
rypbomn = [23.3 21.3 21.0]
rypbosd = [4.0 7.6 7.8]
rynpbo = [11 11 11]

# Rynn Sert
rys1mn = [20.6 12.0 7.8]
rys1sd = [3.6 6.4 5.7]
ryns1 = [11 11 11]

# Rupp pbo - no week 12 here
rpbomn = [17.0 15.3 15.9]
rpbosd = [3.6 4.4 5.3]  # week 12 miscoded?
rnpbo = [63 56 49]

# Rupp flvx
rs1mn = [18.7 11.1 7.1]
rs1sd = [2.9 6.0 6.1]  # week 12 miscoded?
rns1 = [61 52 50]


# Rynn Vfx study pbo - no week 12 here
# N.B: no SDs available beyond week 0
rypbomn = [40.0 30.6 28.53 0.0]
rypbosd = [5.5 5.5 5.5 0.0]
rynpbo = [154 154 154 0]

# Rynn Vfx
rys0mn = [39.75 27.68 24.21 0.0]
rys0sd = [5.0 5.0 5.0 0.0]
ryns0 = [159 159 159 0]


# Strawn Dulox study pbo - week 7 and 10 instead of 8 and 12 here
stpbomn = [17.2 12.59 11.44 11.0]
stpbosd = [2.4 2.4 2.4 2.7]
stnpbo = [133 133 133 133]

# Strawn Dulox
sts0mn = [17.3 10.8 9.6 8.0]
sts0sd = [2.2 2.2 2.2 2.5]
stns0 = [135 135 135 135]


# Geller study pbo ### INCLUDING 12 WEEKS
glpbomn = [16.92 14.09 12.95 13.72]
glpbosd = [1.8 2.3 2.7 2.7]
glnpbo = [80 80 58 58]  # no week 12

#Geller Atom
gls0mn = [17.5 12.33 10.82 12.0]
gls0sd = [1.8 2.3 2.7 2.7]
glns0 = [78 78 55 55]


# March study pbo
mpbomn = [66.2 56.6 51.9 0.0]
mpbosd = [10.6 9.0 12.7 0.0]
mnpbo = [148 148 148 0]  # no week 12

# March Vfx
ms0mn = [64.8 53.1 46.7 0.0]
ms0sd = [10.0 8.9 12.4 0.0]
msns0 = [137 137 137 0]


include("bayes_meta_analysis.jl")
using Plots, StatPlots, PlotRecipes
gr()

##### select study below
# Walkup
wmnp, wsep = standardize_scores(wpbomn,wpbosd,wnpbo)
wmnt, wset = standardize_scores(ws1mn,ws1sd,wns1)

drawst, drawsnt = post_treat_notreat(wmnt, wset.^2, wns1, wmnp, wsep.^2, wnpbo)

# Birmaher
bmnp, bsep = standardize_scores(bpbomn,bpbosd,bnpbo)
bmnt, bset = standardize_scores(bs1mn,bs1sd,bns1)

drawst, drawsnt = post_treat_notreat(bmnt, bset.^2, bns1, bmnp, bsep.^2, bnpbo)

# Wagner
wgmnp, wgsep = standardize_scores(wgpbomn,wgpbosd,wgnpbo)
wgmnt, wgset = standardize_scores(wgs1mn,wgs1sd,wgns1)

drawst, drawsnt = post_treat_notreat(wgmnt, wgset.^2, wgns1, wgmnp, wgsep.^2, wgnpbo)

# Rynn
rymnp, rysep = standardize_scores(rypbomn,rypbosd,rynpbo)
rymnt, ryset = standardize_scores(rys1mn,rys1sd,ryns1)

drawst, drawsnt = post_treat_notreat(rymnt, ryset.^2, ryns1, rymnp, rysep.^2, rynpbo)

# Rupp
rmnp, rsep = standardize_scores(rpbomn,rpbosd,rnpbo)
rmnt, rset = standardize_scores(rs1mn,rs1sd,rns1)

drawst, drawsnt = post_treat_notreat(rmnt, rset.^2, rns1, rmnp, rsep.^2, rnpbo)

# plot study selected
drawst, drawsnt = post_treat_notreat(wmnt, wset.^2, wns1, wmnp, wsep.^2, wnpbo)
plot(drawst[:,1], st=:density, linecolor = "blue" , label="Walkup",gridcolor="grey", fill=(0,:blue),alpha=0.2,title="")
plot!(drawst[:,2], st=:density, linecolor = "blue" , label="",gridcolor="grey", fill=(0,:blue),alpha=0.2,title="")
plot!(drawst[:,3], st=:density, linecolor = "blue" , label="",gridcolor="grey", fill=(0,:blue),alpha=0.2,title="")
#title!("Walkup")


drawst, drawsnt = post_treat_notreat(rmnt, rset.^2, rns1, rmnp, rsep.^2, rnpbo)
plot!(drawst[:,1], st=:density, linecolor = "green" , label="Rupp",gridcolor="grey", fill=(0,:green),alpha=0.2,title="")
plot!(drawst[:,2], st=:density, linecolor = "green" , label="",gridcolor="grey", fill=(0,:green),alpha=0.2,title="")
plot!(drawst[:,3], st=:density, linecolor = "green" , label="",gridcolor="grey", fill=(0,:green),alpha=0.2,title="")
title!("Walkup & Rupp")
savefig("WalkupaandRupp.png")
# Combining Walkup and Birmaher

# placebo update
m1p = wmnp
s21p = wsep.^2
n1p = wnpbo
m0p = bmnp
s20p = bsep.^2
n0p = bnpbo

m2p = zeros(nwks)
s22p = zeros(nwks)
n2p = zeros(nwks)
for i in 1:nwks
  m2p[i], s22p[i], n2p[i] = update_mean(m1p[i],m0p[i],s21p[i],s20p[i],n1p[i],n0p[i])
end

# treatment update
m1t = wmnt
s21t = wset.^2
n1t = wns1
m0t = bmnt
s20t = bset.^2
n0t = bns1

m2t = zeros(nwks)
s22t = zeros(nwks)
n2t = zeros(nwks)
for i in 1:nwks
  m2t[i], s22t[i], n2t[i] = update_mean(m1t[i],m0t[i],s21t[i],s20t[i],n1t[i],n0t[i])
end

# draw from posteriors for treatment and pcbo
drawst, drawsnt = post_treat_notreat(m2t, s22t, n2t, m2p, s22p, n2p)



# Now combining first three studies
# placebo update
m1p = wgmnp
s21p = wgsep.^2
n1p = wgnpbo
m0p = m2p
s20p = s22p
n0p = n2p

m2p = zeros(nwks)
s22p = zeros(nwks)
n2p = zeros(nwks)
for i in 1:nwks
  m2p[i], s22p[i], n2p[i] = update_mean(m1p[i],m0p[i],s21p[i],s20p[i],n1p[i],n0p[i])
end

# treatment update
m1t = wgmnt
s21t = wgset.^2
n1t = wgns1
m0t = m2t
s20t = s22t
n0t = n2t

m2t = zeros(nwks)
s22t = zeros(nwks)
n2t = zeros(nwks)
for i in 1:nwks
  m2t[i], s22t[i], n2t[i] = update_mean(m1t[i],m0t[i],s21t[i],s20t[i],n1t[i],n0t[i])
end

# Now combining first four studies
# placebo update
m1p = rymnp
s21p = rysep.^2
n1p = rynpbo
m0p = m2p
s20p = s22p
n0p = n2p

m2p = zeros(nwks)
s22p = zeros(nwks)
n2p = zeros(nwks)
for i in 1:nwks
  m2p[i], s22p[i], n2p[i] = update_mean(m1p[i],m0p[i],s21p[i],s20p[i],n1p[i],n0p[i])
end

# treatment update
m1t = rymnt
s21t = ryset.^2
n1t = ryns1
m0t = m2t
s20t = s22t
n0t = n2t

m2t = zeros(nwks)
s22t = zeros(nwks)
n2t = zeros(nwks)
for i in 1:nwks
  m2t[i], s22t[i], n2t[i] = update_mean(m1t[i],m0t[i],s21t[i],s20t[i],n1t[i],n0t[i])
end

# Now combining all five  SSRI studies
# placebo update
m1p = rmnp
s21p = rsep.^2
n1p = rnpbo
m0p = m2p
s20p = s22p
n0p = n2p

m2p = zeros(nwks)
s22p = zeros(nwks)
n2p = zeros(nwks)
for i in 1:nwks
  m2p[i], s22p[i], n2p[i] = update_mean(m1p[i],m0p[i],s21p[i],s20p[i],n1p[i],n0p[i])
end

# treatment update
m1t = rmnt
s21t = rset.^2
n1t = rns1
m0t = m2t
s20t = s22t
n0t = n2t

m2t = zeros(nwks)
s22t = zeros(nwks)
n2t = zeros(nwks)
for i in 1:nwks
  m2t[i], s22t[i], n2t[i] = update_mean(m1t[i],m0t[i],s21t[i],s20t[i],n1t[i],n0t[i])
end



# draw from posteriors for treatment and pcbo
drawst, drawsnt = post_treat_notreat(m2t, s22t, n2t, m2p, s22p, n2p)


# Combining Walkup and Wagner

# placebo update
m1p = wmnp
s21p = wsep.^2
n1p = wnpbo
m0p = wgmnp
s20p = wgsep.^2
n0p = wgnpbo

m2p = zeros(nwks)
s22p = zeros(nwks)
n2p = zeros(nwks)
for i in 1:nwks
  m2p[i], s22p[i], n2p[i] = update_mean(m1p[i],m0p[i],s21p[i],s20p[i],n1p[i],n0p[i])
end

# treatment update
m1t = wmnt
s21t = wset.^2
n1t = wns1
m0t = wgmnt
s20t = wgset.^2
n0t = wgns1

m2t = zeros(nwks)
s22t = zeros(nwks)
n2t = zeros(nwks)
for i in 1:nwks
  m2t[i], s22t[i], n2t[i] = update_mean(m1t[i],m0t[i],s21t[i],s20t[i],n1t[i],n0t[i])
end

# draw from posteriors for treatment and pcbo
drawst, drawsnt = post_treat_notreat(m2t, s22t, n2t, m2p, s22p, n2p)


# Plot of treatment group posterior mean densities
plot(drawst[:,1], st=:density, linecolor = 1 , gridcolor=:lightgrey,label="T wk 0")
plot!(drawst[:,2], st=:density, linecolor = 1, label="T wk 4")
plot!(drawst[:,3], st=:density, linecolor = 1, label="T wk 8" )
plot!(drawst[:,4], st=:density, linecolor = 1, label="T wk 12" )

# Plot of pcbo group posterior mean densities
plot!(drawsnt[:,1], st=:density, linecolor = 2 , label="Pbo wk 0")
plot!(drawsnt[:,2], st=:density, linecolor = 2, label="Pbo wk 4")
plot!(drawsnt[:,3], st=:density, linecolor = 2, label="Pbo wk 8" )
plot!(drawsnt[:,4], st=:density, linecolor = 2, label="Pbo wk 12" )

# Plot posterior difference in means
plot((drawst[:,1]-drawsnt[:,1]), st=:density, linecolor = 1 , gridcolor=:lightgrey,label="T wk 0")
plot!((drawst[:,2]-drawsnt[:,2]), st=:density, linecolor = 2, label="T wk 4")
plot!((drawst[:,3]-drawsnt[:,3]), st=:density, linecolor = 3, label="T wk 8" )
plot!((drawst[:,4]-drawsnt[:,4]), st=:density, linecolor = 4, label="T wk 12" )

# now update with the next sample and repeat the plots of the difference in means
# Also do boxplots for each posterior mean difference

# plot mean and 0.95 interval
mns = zeros(nwks)
lub = zeros(nwks,2)
for i in 1:nwks
 mns[i] = mean(drawst[:,i]-drawsnt[:,i])
 lub[i,:] = quantile((drawst[:,i]-drawsnt[:,i]),[0.025 0.975])
end

plot(mns,label="mean",linecolor=1)
plot!(lub[:,1],label="lower 0.95", linecolor=3)
plot!(lub[:,2],label="upper 0.95", linecolor=3)







# how to plot in a loop? Below doesn't work
M, nwks = size(drawst)
plot(drawst[:,1], st=:density, linecolor = 1 , gridcolor=:lightgrey,label=1)
for i in 2:nwks
  plot!(drawst[:,i], st=:density, linecolor = 1, label=i)
end



###############################
##### old stuff below
# Dividing by week 0 score for now:
stdz_wkup_mn = wpbomn/wpbomn[1]
stdz_berm_mn = bpbomn/bpbomn[1]
stdz_rupp_mn = rpbomn/rpbomn[1]

# scaled SD
stdz_wkup_sd = wpbosd/wpbomn[1]
stdz_berm_sd = bpbosd/bpbomn[1]
stdz_rupp_sd = rpbosd/rpbomn[1]

# standard error of scaled mean
stdz_wkup_se = stdz_wkup_sd./sqrt(wnpbo)
stdz_berm_se = stdz_berm_sd./sqrt(bnpbo)
stdz_rupp_se = stdz_rupp_sd./sqrt(rnpbo)

using Distributions, Plots, PlotRecipes
pyplot()


## define function for Bayesian updating of mean
## using Gelman et al. (2014), p.64-68.

### m0, s20, n0 = prior mean, variance, number of obs.
### m1, s21, n1 = sample mean, variance, number of obs.
### If another sample is obtained, m1 becomes m0, etc. to update



mns = [stdz_wkup_mn; stdz_berm_mn; stdz_rupp_mn]
ses = [stdz_wkup_se; stdz_berm_se; stdz_rupp_se]
ns = [wnpbo bnpbo rnpbo]
labels = ["Walkup" "Berm" "Rupp"]


j = 1
i = 1 # week
mpost_sample1 = marg_post_mu(mns[1,1],ses[1,1]^2,ns[1],M)
plot(mpost_sample1, st=:density, linecolor = i , gridcolor=:lightgrey,label=labels[1])

nstudies, nwks = size(mns)
# for i = 1:nstudies
  j = 3
  i = 3
    mpost_sample1 = marg_post_mu(mns[i,j],ses[i,j]^2,ns[i],M)
    plot!(mpost_sample1, st=:density, linecolor = i ,label=labels[i])

# end
