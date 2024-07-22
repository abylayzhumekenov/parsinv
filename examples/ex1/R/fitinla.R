# load libraries
library(INLAspacetime)
library(inlabru)


# command line arguments
args = commandArgs(trailingOnly=FALSE)
for(i in seq_along(args)){
    if(args[i] == "-nt") nt = as.integer(args[i+1])
    if(args[i] == "-ns") ns = as.integer(args[i+1])
}
if(!exists(deparse(substitute(ns)))) ns = 12
if(!exists(deparse(substitute(nt)))) nt = 2
source("R/simulation/generate.R")


# ------------------------------------------------------------------------------


# define a model using INLAspacetime
data = list(xcoord = rep(sloc[,1], nt), ycoord = rep(sloc[,2], nt), zcoord = rep(sloc[,3], nt),
            time = rep(1:nt, each=ms), x = x, y = y)
model = y ~ -1 + Intercept(1) + x + field(list(space = cbind(xcoord, ycoord, zcoord), time = time), model = model.st)
model.st = stModel.define(smesh, tmesh, "121",
                          control.priors = list(prs    = c(1.00, 0.01),
                                                prt    = c(1.00, 0.01),
                                                psigma = c(1.00, 0.01)))
lkprec = list(prec = list(initial = 1.00, fixed = FALSE, prior = "pc.prec", param = c(1.00, 0.01)))


# fit using INLA
result = bru(model, 
             like(formula = y ~ ., 
                  family = "gaussian",
                  control.family = list(hyper = lkprec), 
                  data = data),
             options = list(verbose = TRUE,
                            safe = FALSE,
                            control.inla = list(int.strategy = "eb"),
                            control.fixed = list(prec = list(prec = 1e-5, prec.intercept = 1e-5))))
print(unname(result$misc$configs$config[[1]]$theta[c(2:4,1)]))