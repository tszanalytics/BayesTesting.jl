
# data from Dosing Extracted Data 2-12-16
#
include("load_data_standardize.jl")

## set to 4 to look at studies with 12 week data
# set to 3 for only up to 8 weeks
nwks = 4

# Walkup & Wagner combined:
# Placebo
m2wp = zeros(nwks,1)
s2wp = zeros(nwks,1)
n2wp = zeros(nwks,1)
for i in 1:nwks
  m2wp[i], s2wp[i], n2wp[i] = update_mean(wmnp[i],wgmnp[i],wsep[i]^2,wgsep[i]^2,wnpbo[i],wgnpbo[i])
end

# treatment
m2wt = zeros(nwks,1)
s2wt = zeros(nwks,1)
n2wt = zeros(nwks,1)
for i in 1:nwks
  m2wt[i], s2wt[i], n2wt[i] = update_mean(wmnt[i],wgmnt[i],wset[i]^2,wgset[i]^2,wns1[i],wgns1[i])
end


## Add in other studies
# pcbo
for i in 1:nwks
  # Rynn
  m2wp[i], s2wp[i], n2wp[i] = update_mean(rymnp[i],m2wp[i],rysep[i]^2,s2wp[i],rynpb1[i],n2wp[i])
  # Rupp
  m2wp[i], s2wp[i], n2wp[i] = update_mean(rmnp[i],m2wp[i],rsep[i]^2,s2wp[i],rnpbo[i],n2wp[i])
  # Birmaher
  m2wp[i], s2wp[i], n2wp[i] = update_mean(bmnp[i],m2wp[i],bsep[i]^2,s2wp[i],bnpbo[i],n2wp[i])

end

# treatment
for i in 1:nwks
  #Rynn
  m2wt[i], s2wt[i], n2wt[i] = update_mean(rymnt[i],m2wt[i],ryset[i]^2,s2wt[i],ryns1[i],n2wt[i])
  #Rupp
  m2wt[i], s2wt[i], n2wt[i] = update_mean(rmnt[i],m2wt[i],rset[i]^2,s2wt[i],rns1[i],n2wt[i])
  # Birmaher
  m2wt[i], s2wt[i], n2wt[i] = update_mean(bmnt[i],m2wt[i],bset[i]^2,s2wt[i],bns1[i],n2wt[i])

end

M = 100000 # no. of draws

dfwk1 = zeros(M,nwks)  # SSRIs
for i in 1:nwks
  drawst, drawsnt = post_treat_notreat(m2wt[i], s2wt[i], n2wt[i], m2wp[i], s2wp[i], n2wp[i],M=M)
  dfwk1[:,i] = treat_pcbo_difference(drawst,drawsnt)
end

# plot difference (treat - pcbo)  each week
#plot(dfwk1[:,1], st=:density, linecolor = "black" , label="Dwk0",gridcolor="grey",title="SSRIs, T-Pbo",legend=:topleft)
#plot!(dfwk1[:,2], st=:density, linecolor = "green" , label="Dwk4")
#plot!(dfwk1[:,3], st=:density, linecolor = "purple" , label="Dwk8")
#plot!(dfwk1[:,4], st=:density, linecolor = "orange" , label="Dwk12")

plot(dfwk1[:,1], st=:density, linecolor = "green" , label="SSRI wk0",gridcolor="grey", fill=(0,:green),alpha=0.1,title="SSRIs, T-Pbo",legend=:topleft)
plot!(dfwk1[:,2], st=:density, linecolor = "green" , label="SSRI wk4", fill=(0,:green),alpha=0.3)
plot!(dfwk1[:,3], st=:density, linecolor = "green" , label="SSRI wk8", fill=(0,:green),alpha=0.5)
plot!(dfwk1[:,4], st=:density, linecolor = "green" , label="SSRI wk12", fill=(0,:green),alpha=0.7)



savefig("SSRIdiffsAll.png")

bayespval(dfwk1[:,3])
bayespval(dfwk1[:,4])
mcodds(dfwk1[:,3])
mcodds(dfwk1[:,4])

m = zeros(nwks)
for (i=1:nwks)
  m[i] = mean(dfwk1[:,i])
end
traj_plot(m,nwks=3,overlay=true)



#### SNRIs

# pcbo
m3wp = zeros(nwks,1)
s3wp = zeros(nwks,1)
n3wp = zeros(nwks,1)

for i in 1:nwks
  # Rynn + Strawn
  #m3wp[i], s3wp[i], n3wp[i] = update_mean(ry0mnp[i],stmnp[i],ry0sep[i]^2,stsep[i],rynpb0[i],stnpbo[i])
  # Geller
  #m3wp[i], s3wp[i], n3wp[i] = update_mean(glmnp[i],m3wp[i],glsep[i]^2,s3wp[i],glnpbo[i],n3wp[i])
  # March
#  m3wp[i], s3wp[i], n3wp[i] = update_mean(mmnp[i],m3wp[i],msep[i]^2,s3wp[i],mnpbo[i],n3wp[i])

# Geller + Strawn
  m3wp[i], s3wp[i], n3wp[i] = update_mean(glmnp[i],stmnp[i],glsep[i]^2,stsep[i],glnpbo[i],stnpbo[i])

end

# treatment
m3wt = zeros(nwks,1)
s3wt = zeros(nwks,1)
n3wt = zeros(nwks,1)

for i in 1:nwks
  # Rynn + Strawn
#  m3wt[i], s3wt[i], n3wt[i] = update_mean(ry0mnt[i],stmnt[i],ry0set[i]^2,stset[i],ryns0[i],stns0[i])
  # Geller
#  m3wt[i], s3wt[i], n3wt[i] = update_mean(glmnt[i],m3wt[i],glset[i]^2,s3wt[i],glns0[i],n3wt[i])
  # March
#  m3wt[i], s3wt[i], n3wt[i] = update_mean(mmnt[i],m3wt[i],mset[i]^2,s3wt[i],msns0[i],n3wt[i])

  # Geller + Strawn
  m3wt[i], s3wt[i], n3wt[i] = update_mean(glmnt[i],stmnt[i],glset[i]^2,stset[i],glns0[i],stns0[i])

end

M = 100000 # no. of draws

dfwk0 = zeros(M,nwks)  ### SNRIs
for i in 1:nwks
  drawst, drawsnt = post_treat_notreat(m3wt[i], s3wt[i], n3wt[i], m3wp[i], s3wp[i], n3wp[i],M=M)
  dfwk0[:,i] = treat_pcbo_difference(drawst,drawsnt)
end

nwks=4
msnri = zeros(nwks)
for i=1:nwks
  msnri[i] = mean(dfwk0[:,i])
end
traj_plot(msnri,nwks=3)

# plot difference (treat - pcbo)  each week
plot(dfwk0[:,1], st=:density, linecolor = "blue" , label="SNRI wk0",gridcolor="grey",xlims=(-0.8,0.5), fill=(0,:blue),alpha=0.2,title="SSRIs & SNRIs")
plot!(dfwk0[:,2], st=:density, linecolor = "blue" , label="SNRI wk4", fill=(0,:blue),alpha=0.3)
plot!(dfwk0[:,3], st=:density, linecolor = "blue" , label="SNRI wk8", fill=(0,:blue),alpha=0.4)
plot!(dfwk0[:,4], st=:density, linecolor = "blue" , label="SNRI wk12", fill=(0,:blue),alpha=0.4)

plot!(dfwk1[:,1], st=:density, linecolor = "green" , label="SSRI wk0",gridcolor="grey", fill=(0,:green),alpha=0.2)
plot!(dfwk1[:,2], st=:density, linecolor = "green" , label="SSRI wk4", fill=(0,:green),alpha=0.3)
plot!(dfwk1[:,3], st=:density, linecolor = "green" , label="SSRI wk8", fill=(0,:green),alpha=0.4)
plot!(dfwk1[:,4], st=:density, linecolor = "green" , label="SSRI wk12", fill=(0,:green),alpha=0.4)


savefig("SNRIdiffs12wk.png")

## SSRI and SNRI on same fig

plot(dfwk1[:,2], st=:density, linecolor = "purple" , label="SSRI wk4",title="SSRIs & SNRIs")
plot!(dfwk1[:,3], st=:density, linecolor = "blue" , label="SSRI wk8")
plot!(dfwk1[:,4], st=:density, linecolor = "black" , label="SSRI wk12")

plot!(dfwk0[:,2], st=:density, linecolor = "red" , label="SNRI wk4")
plot!(dfwk0[:,3], st=:density, linecolor = "pink" , label="SNRI wk8")
plot!(dfwk0[:,4], st=:density, linecolor = "orange" , label="SNRI wk12")

savefig("SNRIvSSRI12wk.pgn")

# create trajectory
md1 = mean(dfwk1[:,1])
ci95d1 = quantile(dfwk1[:,1],[0.025,0.975])
md2 = mean(dfwk1[:,2])
ci95d2 = quantile(dfwk1[:,2],[0.025,0.975])
md3 = mean(dfwk1[:,3])
ci95d3 = quantile(dfwk1[:,3],[0.025,0.975])
md4 = mean(dfwk1[:,4])
ci95d4 = quantile(dfwk1[:,4],[0.025,0.975])

mn = [md1 md2 md3 md4]
cilo = [ci95d1[1] ci95d2[1] ci95d3[1] ci95d4[1] ]
ciup = [ci95d1[2] ci95d2[2] ci95d3[2] ci95d4[2]]

wk = [0 4 8 12]

plot(wk[:],cilo[:],st = :line, color=:green, label="0.95 CI SSRI",title="SSRIs")
plot!(wk[:],mn[:],st = :line, color=:black, label="SSRI mean")
plot!(wk[:],ciup[:],st = :line, color=:green, label="")

scatter!(wk[:],mn[:],color=:black, label="")
scatter!(wk[:],cilo[:],color=:green, label="")
scatter!(wk[:],ciup[:],color=:green, label="")
hline!([0.0],color=:red,label="")
#title!("SSRI - SNRI Efficacy", fontsize = 10pt)


md1 = mean(dfwk0[:,1])
ci95d1 = quantile(dfwk0[:,1],[0.025,0.975])
md2 = mean(dfwk0[:,2])
ci95d2 = quantile(dfwk0[:,2],[0.025,0.975])
md3 = mean(dfwk0[:,3])
ci95d3 = quantile(dfwk0[:,3],[0.025,0.975])
md4 = mean(dfwk0[:,4])
ci95d4 = quantile(dfwk0[:,4],[0.025,0.975])


mn = [md1 md2 md3 md4]
cilo = [ci95d1[1] ci95d2[1] ci95d3[1] ci95d4[1] ]
ciup = [ci95d1[2] ci95d2[2] ci95d3[2] ci95d4[2] ]

plot!(wk[:],cilo[:],st = :line, color=:purple, label="0.95 CI SNRI")
plot!(wk[:],mn[:],st = :line, color=:blue, label="SNRI mean")
plot!(wk[:],ciup[:],st = :line, color=:purple, label="")

scatter!(wk[:],mn[:],color=:black, label="")
scatter!(wk[:],cilo[:],color=:green, label="")
scatter!(wk[:],ciup[:],color=:green, label="")
hline!([0.0],color=:red,label="")
title!("SSRI & SNRI Efficacy", fontsize = 10pt)

savefig("SNRIvSSRItrajectory.png")

# how to plot log function for this?
#X = [ones(3,1) wk[:]]
#y = -log(abs(mn[:]))
#b = X \ y
#w = 0:0.1:8.0
#plot(z->b[1]+b[2]*w)

# Density of difference between SSRI and SNRI


dfSRNR = zeros(M,nwks)  ### SNRIs
for i in 1:nwks
  dfSRNR[:,i] = dfwk1[:,i] - dfwk0[:,i]
end
plot(dfSRNR[:,1], st=:density, linecolor = "green" , label="DEffic wk0")
plot!(dfSRNR[:,2], st=:density, linecolor = "purple" , label="DEffic wk4")
plot!(dfSRNR[:,3], st=:density, linecolor = "red" , label="DEffic wk8")
plot!(dfSRNR[:,4], st=:density, linecolor = "orange" , label="DEffic wk12")

title!("Mean difference in SSRI efficacy")
xaxis!("DPARS relative to placebo")
vline!([0.0],color=:black,label="")
savefig("diffeff.png")

## neither GR nor PGF plots seem to work ?!
# try PGFPlots
#pgfplots()
#gr()
d1 = dfSRNR[:,1]
d2 = dfSRNR[:,2]
d3 = dfSRNR[:,3]
d4 = dfSRNR[:,4]


plot(d1,  st=:density, linecolor = "green",label="week 0", lw=2, fill=(0,:green),alpha=0.3)
#ax[:legend](loc="upper center")
plot!(d2,  st=:density, linecolor = "purple", label="week 4", lw=2, fill=(0,:red),alpha=0.3)
plot!(d3,  st=:density, linecolor = "blue", label="week 8", lw=2, fill=(0,:blue),alpha=0.3)
plot!(d4,  st=:density, linecolor = "red", label="week 12", lw=2, fill=(0,:blue),alpha=0.3)
title!("Efficacy difference: SSRI - SNRI")
#savefig("effic-diff3.png")

# create time plot
md1 = mean(d1)
ci95d1 = quantile(d1,[0.025,0.975])
md2 = mean(d2)
ci95d2 = quantile(d2,[0.025,0.975])
md3 = mean(d3)
ci95d3 = quantile(d3,[0.025,0.975])

mn = [md1 md2 md3]
cilo = [ci95d1[1] ci95d2[1] ci95d3[1] ]
ciup = [ci95d1[2] ci95d2[2] ci95d3[2] ]

wk = [0 4 8]

plot(wk[:],cilo[:],st = :line, color=:green, label="0.975")
plot!(wk[:],mn[:],st = :line, color=:black, label="mean")
plot!(wk[:],ciup[:],st = :line, color=:green, label="0.025")

scatter!(wk[:],mn[:],color=:black, label="")
scatter!(wk[:],cilo[:],color=:green, label="")
scatter!(wk[:],ciup[:],color=:green, label="")
hline!([0.0],color=:red,label="")
title!("SSRI - SNRI Efficacy", fontsize = 10pt)
savefig("SSRI-SNRI3.pgn")


scatter!(u,p)
# not working
#pgfplots()
plot(d1,  st=:density, linecolor = "green",label="week 0", lw=3, fill=(0,:green),alpha=0.3)
plot!(d2,  st=:density, linecolor = "red", label="week 4", lw=3, fill=(0,:red),alpha=0.3)
plot!(d3,  st=:density, linecolor = "blue", label="week 8", lw=3, fill=(0,:blue),alpha=0.3)
plot!(d4,  st=:density, linecolor = "yellow", label="week 12", lw=3, fill=(0,:yellow),alpha=0.3)
vline!([0.0],color="black",label="")
title!("Efficacy difference: SSRI - SNRI")
savefig("effic-diff4.png")
