# LOAD DATA, standardize scores, load Julia packages
# data from Dosing Extracted Data 2-12-16
# Walkup pbo
week = [0 4 8 12]
# nwks = length(week)
### use only up to week 8
nwks = 4

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

# Rynn pbo w Sert study- no week 12 here
rypb1mn = [23.3 21.3 21.0 0.0]
rypb1sd = [4.0 7.6 7.8 0.0]
rynpb1 = [11 11 11 0]

# Rynn Sert
rys1mn = [20.6 12.0 7.8 0.0]
rys1sd = [3.6 6.4 5.7 0.0]
ryns1 = [11 11 11 0]

#### Rynn vfx-pcbo different measure
#rypb2mn = [40.0 30.60 28.51 0.0]
#rypb2sd = [4.0 4.0 4.0 0.0]
#rynpb2 = [159 159 159 0]


# Rupp pbo - no week 12 here
rpbomn = [17.0 15.3 15.9 0.0]
rpbosd = [3.6 4.4 5.3 0.0]  # week 12 miscoded?
rnpbo = [63 56 49 0]

# Rupp flvx
rs1mn = [18.7 11.1 7.1 0.0]
rs1sd = [2.9 6.0 6.1 0.0]  # week 12 miscoded?
rns1 = [61 52 50 0]


# Rynn Vfx study pbo - no week 12 here
# N.B.: no SDs available beyond week 0
rypb0mn = [40.0 30.6 28.53 0.0]
rypb0sd = [4.3 4.3 4.3 0.0]
rynpb0 = [154 154 154 0]

# Rynn Vfx
rys0mn = [39.75 27.68 24.21 0.0]
rys0sd = [3.8 3.8 3.8 0.0]
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
include("CSV_file_creation.jl")

using Plots, StatPlots, PlotRecipes
gr()

# High dose SSRIs

# Walkup
wmnp, wsep = standardize_scores(wpbomn,wpbosd,wnpbo)
wmnt, wset = standardize_scores(ws1mn,ws1sd,wns1)

# Wagner
wgmnp, wgsep = standardize_scores(wgpbomn,wgpbosd,wgnpbo)
wgmnt, wgset = standardize_scores(wgs1mn,wgs1sd,wgns1)

# Rupp
rmnp, rsep = standardize_scores(rpbomn,rpbosd,rnpbo)
rmnt, rset = standardize_scores(rs1mn,rs1sd,rns1)

# low dose SSRIs: Birhamer Rynn

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
