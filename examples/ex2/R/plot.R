# load the libraries
suppressMessages(suppressWarnings(library(INLA)))
library(ggplot2)
library(maps)
library(mapdata)
library(sp)
library(INLA)
source("parsinv.petsc.io.R")

# load the saved data
load("data/smesh.Rdata")
load("data/inla.Rdata")
theta.parsinv = parsinv.vec.read("data/theta")
muu.parsinv = parsinv.vec.read("data/muu")
mub.parsinv = parsinv.vec.read("data/mub")
sdu.parsinv = parsinv.vec.read("data/sdu")
sdb.parsinv = parsinv.vec.read("data/sdb")



# create a US border
usa = map_data("usa")
usa = usa[1:6886,]
coordinates(usa) = ~ long + lat
usa@proj4string = CRS("+proj=longlat +datum=WGS84 +lon_0=100")
usa = spTransform(usa, CRS("+proj=moll +units=km"))
borderline = usa@coords
borderline[,1] = (borderline[,1]-min(bound$loc[,1])) / (diff(range(bound$loc[,1])))
borderline[,2] = (borderline[,2]-min(bound$loc[,2])) / (diff(range(bound$loc[,2])))

# create a projection
grid.ratio = diff(range(bound$loc[,2])) / diff(range(bound$loc[,1]))
grid.res = c(256, 256)
grid.res[2] = round(grid.res[1] * grid.ratio)
grid.x = seq(min(bound$loc[,1]), max(bound$loc[,1]), length = grid.res[1])
grid.y = seq(min(bound$loc[,2]), max(bound$loc[,2]), length = grid.res[2])
grid.loc = expand.grid(grid.x, grid.y)
grid.loc = as.matrix(grid.loc)
grid.A = inla.spde.make.A(mesh = smesh, loc = grid.loc)

# set the plot range
zlim.muu = range(muu.parsinv[1:(ns*nt)])
zlim.sdu = range(sdu.parsinv[1:(ns*nt)])



# plot several slices
tt = c(21, 22, 365, 4383)
pdf("data/fig.app.1.pdf", width=14, height=5)
par(mfcol=c(2,length(tt)), mar=c(2,2,0,0))
for(i in seq_along(tt)){
  # set time
  t = tt[i]
  
  # project the mean
  muu.parsinv.grid = drop(grid.A %*% muu.parsinv[1:ns + (t-1)*ns])
  muu.parsinv.grid = matrix(muu.parsinv.grid, grid.res[1])
  
  # project the sd
  sdu.parsinv.grid = drop(grid.A %*% sdu.parsinv[1:ns + (t-1)*ns])
  sdu.parsinv.grid = matrix(sdu.parsinv.grid, grid.res[1])
  
  # plot the spatial field
  image(muu.parsinv.grid, col=viridisLite::viridis(100), asp=grid.ratio, zlim=zlim.muu, xaxt="n", yaxt="n", axes=FALSE)
  lines(borderline, col="white", lwd=2)
  if(i==1) title(ylab="Mean", line=1)
  image(sdu.parsinv.grid, col=viridisLite::inferno(100), asp=grid.ratio, zlim=zlim.sdu, xaxt="n", yaxt="n", axes=FALSE)
  lines(borderline, col="white", lwd=2)
  if(i==1) title(ylab="SD", line=1)
  title(xlab=paste0("t = ", t), line = 1)
}
dev.off()

n_col = 100
pdf("data/fig.app.1.colorbar.pdf", width=1, height=5)
par(mfcol=c(2,1), mar=c(2,2.5,1,1))
image(x = 1, 
      y = seq(zlim.muu[1], zlim.muu[2], length=n_col), 
      z = matrix(seq(zlim.muu[1], zlim.muu[2], length=n_col), nrow=1),
      col = viridisLite::viridis(n_col), axes=FALSE, xlab=NA, ylab=NA)
axis(side=2, tick=FALSE, at=zlim.muu, line=-1, labels=round(zlim.muu, 1), las=1)
image(x = 1, 
      y = seq(zlim.sdu[1], zlim.sdu[2], length=n_col), 
      z = matrix(seq(zlim.sdu[1], zlim.sdu[2], length=n_col), nrow=1), 
      col = viridisLite::inferno(n_col), axes=FALSE, xlab=NA, ylab=NA)
axis(side=2, tick=FALSE, at=zlim.sdu, line=-1, labels=round(zlim.sdu, 1), las=1)
dev.off()