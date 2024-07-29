# load libraries and helper functions
library(INLA)
library(fmesher)
source("parsinv.stmatrices.R")
source("parsinv.stmodel.R")
source("parsinv.petsc.io.R")


# command line arguments
args = commandArgs(trailingOnly=FALSE)
for(i in seq_along(args)){
    if(args[i] == "-nt") nt = as.integer(args[i+1])
    if(args[i] == "-ns") ns = as.integer(args[i+1])
}
if(!exists(deparse(substitute(ns)))) ns = 12
if(!exists(deparse(substitute(nt)))) nt = 2
res = round(sqrt((ns-2)/10))


# ------------------------------------------------------------------------------


# create meshes
tmesh = fm_mesh_1d(1:nt)
smesh = fm_rcdt_2d_inla(globe = res)
nt = tmesh$n
ns = smesh$n


# simulate locations
set.seed(1)
ms = 10*ns
sloc = matrix(rnorm(ms*3), ms)
sloc = sloc / sqrt(rowSums(sloc^2))
amat = parsinv.amatrices(tmesh, smesh, NULL, sloc)
A = t(kronecker(amat$At, amat$As))


# simulate (unstructured) data
u = sin(rep(1:nt, each = ns)/nt*2*pi*3) + rnorm(nt*ns)
x = rep(1:nt, each = ms)
y = 1 + x + drop(A%*%u) + rnorm(length(x), 0, 0.1)


# ------------------------------------------------------------------------------


# save the model information
stmodel = parsinv.stmodel.define("121", tmesh$manifold, smesh$manifold)
parsinv.stmodel.write(stmodel, "../data/stmodel")


# save the data
Ab = cbind(1, x)
parsinv.dense.write(Ab, "../data/Ab")
parsinv.vec.write(y, "../data/y")


# save stmatrices
jmat = parsinv.jmatrices(stmodel, tmesh)
gmat = parsinv.gmatrices(stmodel, smesh)
# amat = parsinv.amatrices(tmesh, smesh)
parsinv.mats.write(jmat, "../data")
parsinv.mats.write(gmat, "../data")
parsinv.mats.write(amat, "../data")