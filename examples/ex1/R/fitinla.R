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
if(!exists(deparse(substitute(nt)))) nt = 4
source("generate.R")


# ------------------------------------------------------------------------------


# define a model using INLAspacetime
data = list(xcoord = rep(sloc[,1], nt), ycoord = rep(sloc[,2], nt), zcoord = rep(sloc[,3], nt),
            time = rep(1:nt, each=ms), x = x, y = y)
model = y ~ -1 + Intercept(1) + x + field(list(space = cbind(xcoord, ycoord, zcoord), time = time), model = model.st)
model.st = stModel.define(smesh, tmesh, "121",
                          control.priors = list(prs    = c(1.00, 0.00),
                                                prt    = c(10.0, 0.00),
                                                psigma = c(1.00, 0.00)))  
                          # in original scale; change to c(1.00, 0.01) to optimize theta
lkprec = list(prec = list(initial = 1.00, fixed = TRUE, prior = "pc.prec", param = c(1.00, 0.01)))
                          # in log scale; change to fixed = FALSE to optimize theta

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

# save results
theta.inla = unname(result$misc$configs$config[[1]]$theta[c(2:4,1)])  # empty, when theta is fixed!
muu.inla = result$summary.random$field$mean
mub.inla = result$summary.fixed$mean
sdu.inla = result$summary.random$field$sd
sdb.inla = result$summary.fixed$sd
save(list = c("theta.inla", "muu.inla", "mub.inla", "sdu.inla", "sdb.inla"), file="data/inla.Rdata")
