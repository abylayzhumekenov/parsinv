# load libraries and helper functions
library(INLA)
library(fmesher)
library(INLAspacetime)
library(inlabru)
source("R/parsinv.stmatrices.R")
source("R/parsinv.stmodel.R")
source("R/parsinv.petsc.io.R")


# command line arguments
args = commandArgs(trailingOnly=FALSE)
for(i in seq_along(args)){
    if(args[i] == "-nt") nt = as.integer(args[i+1])
    if(args[i] == "-ns") ns = as.integer(args[i+1])
    if(args[i] == "-hh") theta = as.double(args[i+1:4])
}
if(!exists(deparse(substitute(ns)))) ns = 12
if(!exists(deparse(substitute(nt)))) nt = 2
if(!exists(deparse(substitute(theta)))) theta = rep(0, 4)
res = round(sqrt((ns-2)/10))


# ------------------------------------------------------------------------------


# create meshes
tmesh = fm_mesh_1d(1:nt)
smesh = fm_rcdt_2d_inla(globe = res)
nt = tmesh$n
ns = smesh$n


# simulate locations
set.seed(1)
ms = 3*ns
sloc = matrix(rnorm(ms*3), ms)
sloc = sloc / sqrt(rowSums(sloc^2))
amat = parsinv.amatrices(tmesh, smesh, NULL, sloc)
A = t(kronecker(amat$At, amat$As))


# ------------------------------------------------------------------------------


# simulate from the SPDE using INLA
data = data.frame(xcoord = 1, ycoord = 0, zcoord = 0, time = 1, x = 0, y = NA)
model = y ~ -1 + Intercept(1) + x + field(list(space = cbind(xcoord, ycoord, zcoord), time = time), model = model.st)
model.st = stModel.define(smesh, tmesh, "121",
                          control.priors = list(prs    = c(1.00, 0.01),
                                                prt    = c(1.00, 0.01),
                                                psigma = c(1.00, 0.01)))
lkprec = list(prec = list(initial = 1.00, fixed = FALSE, prior = "pc.prec", param = c(1.00, 0.01)))
int.design = t(c(theta[c(4,1:3)], 1))
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = TRUE,
                            safe = FALSE,
                            control.inla = list(int.strategy = "user", int.design = int.design, mode.known = TRUE),
                            control.fixed = list(prec = list(prec = 1e-5, prec.intercept = 1e-5))))


# simulate the data
u = head(inla.qsample(1, result$misc$config$config[[1]]$Q), -2)
x = rep(1:nt, each = ms)
y = 1 + x + drop(A%*%u) + rnorm(length(x), 0, exp(-theta[4]/2))


# ------------------------------------------------------------------------------


# save the model information
stmodel = parsinv.stmodel.define("121", tmesh$manifold, smesh$manifold)
parsinv.stmodel.write(stmodel, "data/stmodel")


# save the data
Ab = cbind(1, x)
parsinv.dense.write(Ab, "data/Ab")
parsinv.vec.write(y, "data/y")


# save stmatrices
jmat = parsinv.jmatrices(stmodel, tmesh)
gmat = parsinv.gmatrices(stmodel, smesh)
# amat = parsinv.amatrices(tmesh, smesh)
parsinv.mats.write(jmat, "data")
parsinv.mats.write(gmat, "data")
parsinv.mats.write(amat, "data")

print(result$misc$configs$config[[1]]$theta[c(2:4,1)])
