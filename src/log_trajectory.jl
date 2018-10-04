# data from Dosing Extracted Data 2-12-16
# Walkup pbo
week = [0 4 8 12]
# nwks = length(week)

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
rypb1mn = [23.3 21.3 21.0]
rypb1sd = [4.0 7.6 7.8]
rynpb1 = [11 11 11]

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
# N.B.: no SDs available beyond week 0
rypb0mn = [40.0 30.6 28.53 0.0]
rypb0sd = [5.5 5.5 5.5 0.0]
rynpb0 = [154 154 154 0]

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
include("bayes_odds_pvalue.jl")
using Plots, StatPlots, PlotRecipes
gr()

# Walkup
wmnp, wsep = standardize_scores(wpbomn,wpbosd,wnpbo)
wmnt, wset = standardize_scores(ws1mn,ws1sd,wns1)

# Wagner
wgmnp, wgsep = standardize_scores(wgpbomn,wgpbosd,wgnpbo)
wgmnt, wgset = standardize_scores(wgs1mn,wgs1sd,wgns1)

# Rupp
rmnp, rsep = standardize_scores(rpbomn,rpbosd,rnpbo)
rmnt, rset = standardize_scores(rs1mn,rs1sd,rns1)

# low dose combined: Birhamer Rynn

# Rynn
rymnp, rysep = standardize_scores(rypb1mn,rypb1sd,rynpb1)
rymnt, ryset = standardize_scores(rys1mn,rys1sd,ryns1)

# Birmaher
bmnp, bsep = standardize_scores(bpbomn,bpbosd,bnpbo)
bmnt, bset = standardize_scores(bs1mn,bs1sd,bns1)

#### SNRIs

# Rynn
ry0mnp, ry0sep = standardize_scores(rypb0mn,rypb0sd,rynpb0)
ry0mnt, ry0set = standardize_scores(rys0mn,rys0sd,ryns0)

# Strawn
stmnp, stsep = standardize_scores(stpbomn,stpbosd,stnpbo)
stmnt, stset = standardize_scores(sts0mn,sts0sd,stns0)

# Geller
glmnp, glsep = standardize_scores(glpbomn,glpbosd,glnpbo)
glmnt, glset = standardize_scores(gls0mn,gls0sd,glns0)

# March
mmnp, msep = standardize_scores(mpbomn,mpbosd,mnpbo)
mmnt, mset = standardize_scores(ms0mn,ms0sd,msns0)


### nwks = 3 => use only up to week 8
### nwks = 4 => use up to week 12
nwks = 3
# log trajectory model
### want diff(treat - pcbo) means here:
mn = wmnp  # means for trajectory plot
z = [1 2 3]'   # for up to week 8 = 3 data points
# z = [1 2 3 4]'   # for up to week 12 = 4 data points
lz = log.(z)
Xloglin = hcat(ones(nwks),lz)
b = Xloglin \ mn[1:nwks]
s2 = sum((mn[1:nwks] - Xloglin*b).^2)/2
println("coeffs = ",b)
println("s^2 (eqn. variance) = ",s2)

# compute R^2
tss = sum((mn[1:nwks]' - mean(mn[1:nwks])).^2)
R2 = 1 - s2*2/tss
println("Rsquared = ",R2)
zs = 1.0:0.1:4.2
score = b[1] + b[2].*log.(zs')

plot(1.0:0.1:4.2,score',size=(400,200),label="")
mns = mn[1:nwks]
scatter!(z',mns',label="")

s = sqrt(s2)
tcrit5 = 4.3  # 5% crit. value from Student-t with 2 d.f.
up95 = score .+ s*tcrit5
lo95 = score .- s*tcrit5

plot(1.0:0.1:4.2,score',size=(400,200),label="")
scatter!(z',mns',label="")
plot!(1.0:0.1:4.2,up95',size=(400,200),label="",color=:green)
plot!(1.0:0.1:4.2,lo95',size=(400,200),label="",color=:green)
hline!([1.0], color=:red,label="")

# allow nwks to vary (add a loop)
# change x axis labels to weeks instead of period
function traj_plot(means;overlay = false,nwks=3)
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

  if (overlay == false)
    plot(1.0:0.1:4.2,score',size=(400,200),label="")
  else
    plot!(1.0:0.1:4.2,score',size=(400,200),label="")
  end
  scatter!(z',mns',label="")
  plot!(1.0:0.1:4.2,up95',size=(400,200),label="",color=:green)
  plot!(1.0:0.1:4.2,lo95',size=(400,200),label="",color=:green)
#  hline!([1.0], color=:red,label="")
end

nwks=4
traj_plot(wmnp,nwks=nwks)
traj_plot(wmnt, overlay=true,nwks=nwks)

# dufference between pcbo and treatment

# Walkup
drawst, drawsnt = post_treat_notreat(wmnt, wset.^2, wns1, wmnp, wsep.^2, wnpbo,M=10000)
dtreatpcbo = treat_pcbo_difference(drawst,drawsnt)

nwks = 4

mndf = zeros(nwks)
for i = 1:nwks
  mndf[i] = mean(dtreatpcbo[:,i])
end
sddf = std(dtreatpcbo[:,nwks])

traj_plot(mndf,nwks=nwks)



# Birmaher
drawst, drawsnt = post_treat_notreat(bmnt, bset.^2, bns1, bmnp, bsep.^2, bnpbo,M=10000)
dtreatpcbo = treat_pcbo_difference(drawst,drawsnt)

mndf = zeros(nwks)
for i = 1:nwks
  mndf[i] = mean(dtreatpcbo[:,i])
end
sddf = std(dtreatpcbo[:,nwks])

traj_plot(mndf,nwks=nwks)


# Rupp
drawst, drawsnt = post_treat_notreat(rmnt, rset.^2, rns1, rmnp, rsep.^2, rnpbo,M=10000)
dtreatpcbo = treat_pcbo_difference(drawst,drawsnt)

nwks = 3
mndf = zeros(nwks)
for i = 1:nwks
  mndf[i] = mean(dtreatpcbo[:,i])
end
sddf = std(dtreatpcbo[:,nwks])

traj_plot(mndf)



##############################################################

##  Combine studies, week by week - Walkup + Wagner
#

# Placebo
m2wp = zeros(nwks,1)
s2wp = zeros(nwks,1)
n2wp = zeros(nwks,1)
for i in 1:nwks
  m2wp[i], s2wp[i], n2wp[i] = update_mean(wmnp[i],wgmnp[i],wsep[i]^2,wgsep[i]^2,wnpbo[i],wgnpbo[i])
end

# high dose pcbo combined: Walkup and Wagner:
himp = m2wp
his2p = s2wp
hinp = n2wp

# treatment
m2wt = zeros(nwks,1)
s2wt = zeros(nwks,1)
n2wt = zeros(nwks,1)
for i in 1:nwks
  m2wt[i], s2wt[i], n2wt[i] = update_mean(wmnt[i],wgmnt[i],wset[i]^2,wgset[i]^2,wns1[i],wgns1[i])
end

# high dose combined: Walkup and Wagner:
himt = m2wt
his2t = s2wt
hint = n2wt


# low dose combined: Birhamer Rynn

# Rynn
rymnp, rysep = standardize_scores(rypb1mn,rypb1sd,rynpb1)
rymnt, ryset = standardize_scores(rys1mn,rys1sd,ryns1)

# Birmaher
bmnp, bsep = standardize_scores(bpbomn,bpbosd,bnpbo)
bmnt, bset = standardize_scores(bs1mn,bs1sd,bns1)

lowmp = zeros(nwks,1)
los2p = zeros(nwks,1)
lownp = zeros(nwks,1)
lowmt = zeros(nwks,1)
los2t = zeros(nwks,1)
lownt = zeros(nwks,1)
for i in 1:nwks
  # Rynn-Birhamer pcbo
  lowmp[i], los2p[i], lownp[i] = update_mean(rymnp[i],bmnp[i],rysep[i]^2,bsep[i]^2,rynpb1[i],bnpbo[i])
  #Rynn-Birhamer treatment
  lowmt[i], los2t[i], lownt[i] = update_mean(rymnt[i],bmnt[i],ryset[i]^2,bset[i]^2,ryns1[i],bns1[i])
end


M = 100000 # no. of draws
dfwkh = zeros(M,nwks)  # High dose studies
dfwkl = zeros(M,nwks)  # Low dose studies
for i in 1:nwks
  drawst, drawsnt = post_treat_notreat(himt[i], his2t[i], hint[i], himp[i], his2p[i], hinp[i],M=M)
  dfwkh[:,i] = treat_pcbo_difference(drawst,drawsnt)
  drawst, drawsnt = post_treat_notreat(lowmt[i], los2t[i], lownt[i], lowmp[i], los2p[i], lownp[i],M=M)
  dfwkl[:,i] = treat_pcbo_difference(drawst,drawsnt)
end


#using Plots, StatPlots, PlotRecipes
#gr()

plot(dfwkh[:,1], st=:density, linecolor = "blue" , label="hi wk0",gridcolor="grey",xlims=(-0.8,0.5), fill=(0,:blue),alpha=0.2,title="High & Low dosage")
plot!(dfwkh[:,2], st=:density, linecolor = "blue" , label="hi wk4", fill=(0,:blue),alpha=0.3)
plot!(dfwkh[:,3], st=:density, linecolor = "blue" , label="hi wk8", fill=(0,:blue),alpha=0.4)
#plot!(dfwk[4], st=:density, linecolor = "purple" , label="Δwk12")

plot!(dfwkl[:,1], st=:density, linecolor = "green" , label="lo wk0",gridcolor="grey", fill=(0,:green),alpha=0.2)
plot!(dfwkl[:,2], st=:density, linecolor = "green" , label="lo wk4", fill=(0,:green),alpha=0.3)
plot!(dfwkl[:,3], st=:density, linecolor = "green" , label="lo wk8", fill=(0,:green),alpha=0.4)
#plot!(dfwk[4], st=:density, linecolor = "purple" , label="Δwk12")
savefig("hilodosage.png")

# Birmaher results exhibit a decline in pcbo
# effectiveness from week 4 to week 8
# Birmaher alone:
M = 100000 # no. of draws
dfwkb = zeros(M,nwks)  # High dose studies
#dfwkl = zeros(M,nwks)  # Low dose studies
for i in 1:nwks
  drawst, drawsnt = post_treat_notreat(bmnt[i], bset[i].^2, bns1[i], bmnp[i], bsep[i].^2, bnpbo[i],M=M)
  dfwkb[:,i] = treat_pcbo_difference(drawst,drawsnt)
#  drawst, drawsnt = post_treat_notreat(lowmt[i], los2t[i], lownt[i], lowmp[i], los2p[i], lownp[i],M=M)
#  dfwkl[:,i] = treat_pcbo_difference(drawst,drawsnt)
end

plot(dfwkb[:,1], st=:density, linecolor = "black" , label="hi wk0",gridcolor="grey")
plot!(dfwkb[:,2], st=:density, linecolor = "green" , label="hi wk4")
plot!(dfwkb[:,3], st=:density, linecolor = "red" , label="hi wk8")
#pl
