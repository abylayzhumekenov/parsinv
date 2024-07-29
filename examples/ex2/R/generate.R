# command line arguments
args = commandArgs(trailingOnly=TRUE)
for(i in seq_along(args)){
    if(args[i] == "-ns") ns = as.integer(args[i+1])
    if(args[i] == "-nt") nt = as.integer(args[i+1])
    if(args[i] == "-ms") ms = as.integer(args[i+1])
    if(args[i] == "-mt") mt = as.integer(args[i+1])
    if(args[i] == "-res1") res1 = as.double(args[i+1])
    if(args[i] == "-res2") res2 = as.double(args[i+1])
    if(args[i] == "-res3") res3 = as.double(args[i+1])
}

# ------------------------------------------------------------------------------

# load data
suppressMessages(suppressWarnings(library(INLA)))
source("getdata.R")
detach("package:data.table", unload = TRUE)
rm(list = setdiff(ls(), c("wdat", "stations", "ns", "nt", "ms", "mt", "res1", "res2", "res3")))

# load libraries and helper functions
library(INLA)
library(fmesher)
source("parsinv.stmatrices.R")
source("parsinv.stmodel.R")
source("parsinv.petsc.io.R")
set.seed(1)

# set mesh resolution and data dimensions
if(!exists(deparse(substitute(res1)))) res1 = 100
if(!exists(deparse(substitute(res2)))) res2 = 200
if(!exists(deparse(substitute(res3)))) res3 = 1000
if(!exists(deparse(substitute(ms)))) ms = dim(wdat)[1]
if(!exists(deparse(substitute(mt)))) mt = dim(wdat)[2]-1

# ------------------------------------------------------------------------------

# create mesh
bound = inla.nonconvex.hull(points = stations@coords, convex = res1, concave = res1, resolution = 100)
smesh = inla.mesh.2d(boundary = bound, max.edge = c(1, 4)*res2, offset = c(1e-3, res3), cutoff = res2/4)
tmesh = inla.mesh.1d(1:mt)
ns = smesh$n
nt = tmesh$n

# create fem objects
fem.t = inla.mesh.fem(tmesh, order = 2)
fem.s = inla.mesh.fem(smesh, order = 3)

# covariates
wdat = wdat[sort(sample(1:(dim(wdat)[1]), ms)), 1:(mt+1)]
sloc = stations@coords[match(wdat$station, stations$station),]
elev = stations$elevation[match(wdat$station, stations$station)]
Ab = cbind(intercept  = rep(1, ms*mt),
            elevation = rep(elev, mt) / 1000,
            northing  = rep(sloc[,2], mt) / 10000,
            sin       = rep(sin(2*pi*(1:mt-1)/365.25), each = ms),
            cos       = rep(cos(2*pi*(1:mt-1)/365.25), each = ms))
nb = ncol(Ab)
y = c(as.matrix(wdat[,-1])) / 10

# ------------------------------------------------------------------------------

# save the model
stmodel = parsinv.stmodel.define("121", tmesh$manifold, smesh$manifold)
parsinv.stmodel.write(stmodel, "../data/stmodel")

# save the data
parsinv.dense.write(Ab, "../data/Ab")
parsinv.vec.write(y, "../data/y")

# generate and save stmatrices
jmat = parsinv.jmatrices(stmodel, tmesh)
gmat = parsinv.gmatrices(stmodel, smesh)
amat = parsinv.amatrices(tmesh, smesh, NULL, sloc)
parsinv.mats.write(jmat, "../data")
parsinv.mats.write(gmat, "../data")
parsinv.mats.write(amat, "../data")

cat("ms = ", ms, "\n",
    "mt = ", mt, "\n",
    "ns = ", ns, "\n",
    "nt = ", nt, "\n",
    "nb = ", nb, "\n",
    sep="")
