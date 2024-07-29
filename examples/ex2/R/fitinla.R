# load libraries
library(INLAspacetime)
library(inlabru)


# command line arguments
args = commandArgs(trailingOnly=FALSE)
for(i in seq_along(args)){
    if(args[i] == "-nt") nt = as.integer(args[i+1])
    if(args[i] == "-ns") ns = as.integer(args[i+1])
    if(args[i] == "-ms") ms = as.integer(args[i+1])
    if(args[i] == "-mt") mt = as.integer(args[i+1])
    if(args[i] == "-res1") res1 = as.double(args[i+1])
    if(args[i] == "-res2") res2 = as.double(args[i+1])
    if(args[i] == "-res3") res3 = as.double(args[i+1])
}
if(!exists(deparse(substitute(ns)))) ms = 10
if(!exists(deparse(substitute(nt)))) mt = 4
if(!exists(deparse(substitute(res1)))) res1 = 200
if(!exists(deparse(substitute(res2)))) res2 = 500
if(!exists(deparse(substitute(res3)))) res3 = 1000
source("generate.R")


# ------------------------------------------------------------------------------


# define a model using INLAspacetime
data = list(xcoord = rep(sloc[,1], nt), ycoord = rep(sloc[,2], nt), time = rep(1:nt, each=ms), 
            elevation = Ab[,2], northing = Ab[,3], sin = Ab[,4], cos = Ab[,5])
model = y ~ -1 + Intercept(1) + elevation + northing + sin + cos + 
    field(list(space = cbind(xcoord, ycoord), time = time), model = model.st)
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
