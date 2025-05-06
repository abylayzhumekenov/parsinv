# load the saved data
source("parsinv.petsc.io.R")
load("data/smesh.Rdata")
load("data/inla.Rdata")
theta.parsinv = parsinv.vec.read("data/theta")
muu.parsinv = parsinv.vec.read("data/muu")
mub.parsinv = parsinv.vec.read("data/mub")
sdu.parsinv = parsinv.vec.read("data/sdu")
sdb.parsinv = parsinv.vec.read("data/sdb")

# create a projection
plot.res = c(256, 128)
plot.asp = plot.res[2] / plot.res[1]
coords.spherical = expand.grid(seq(0,2*pi,length=plot.res[1]+1)[-1], seq(0,pi,length=plot.res[2]+1)[-1])
coords.cartesian = data.frame(x = sin(coords.spherical$Var2) * cos(coords.spherical$Var1),
                              y = sin(coords.spherical$Var2) * sin(coords.spherical$Var1),
                              z = cos(coords.spherical$Var2))
A.spherical = inla.spde.make.A(smesh, as.matrix(coords.cartesian))

for(t in 1:nt){
  # project the mean
  muu.inla.spherical = drop(A.spherical %*% muu.inla[1:ns+(t-1)*ns])
  muu.inla.spherical = matrix(muu.inla.spherical, plot.res[1])
  muu.parsinv.spherical = drop(A.spherical %*% muu.parsinv[1:ns+(t-1)*ns])
  muu.parsinv.spherical = matrix(muu.parsinv.spherical, plot.res[1])
  
  # project the sd
  sdu.inla.spherical = drop(A.spherical %*% sdu.inla[1:ns+(t-1)*ns])
  sdu.inla.spherical = matrix(sdu.inla.spherical, plot.res[1])
  sdu.parsinv.spherical = drop(A.spherical %*% sdu.parsinv[1:ns+(t-1)*ns])
  sdu.parsinv.spherical = matrix(sdu.parsinv.spherical, plot.res[1])
  
  # set the plot range
  zlim.muu = range(muu.inla.spherical, muu.parsinv.spherical)
  zlim.sdu = range(sdu.inla.spherical, sdu.parsinv.spherical)
  
  # plot the spatial field
  par(mfrow=c(2,2), mar=c(2,2,1,1))
  image(muu.inla.spherical, col=viridisLite::viridis(100), asp=plot.asp, zlim=zlim.muu, xaxt="n", yaxt="n", axes=FALSE)
  title(ylab="Mean", line=1)
  image(muu.parsinv.spherical, col=viridisLite::viridis(100), asp=plot.asp, zlim=zlim.muu, xaxt="n", yaxt="n", axes=FALSE)
  image(sdu.inla.spherical, col=viridisLite::inferno(100), asp=plot.asp, zlim=zlim.sdu, xaxt="n", yaxt="n", axes=FALSE)
  title(xlab="R-INLA", ylab="SD", line=1)
  image(sdu.parsinv.spherical, col=viridisLite::inferno(100), asp=plot.asp, zlim=zlim.sdu, xaxt="n", yaxt="n", axes=FALSE)
  title(xlab="Ovelapping RBMC", line=1)
}