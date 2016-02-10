#  Waterpump parameterisation program.
#  Load 3 sets of data from TOT111, TOT112, and MAT102
#  Input parameters, read in from a data source.


setwd("C:/Users/warwi/OneDrive/R Language/WaterPump3D")

Pumpinfo <- read.csv ("WaterPumpPerfData.csv")
TOT111 <- read.csv ("C:/Users/warwi/Onedrive/Datasets/TOT111 1h.csv")
TOT112 <- read.csv ("C:/Users/warwi/Onedrive/Datasets/TOT112 1h.csv")
MAT102 <- read.csv ("C:/Users/warwi/Onedrive/Datasets/MAT102 1h.csv")

Pumpinfo
PumpNo=5

# read the data for the pump then solve the simultaneous equation for a,b,c where h=a+bq+cq^2  (q units= m^3/s)  
RPM1 <-Pumpinfo[PumpNo,3] # Baseline RPM
beta1 <-Pumpinfo[PumpNo,4]  #h1
beta2 <-Pumpinfo[PumpNo,6]  #h2
beta3 <-Pumpinfo[PumpNo,8]  #h3
alpha1 <-Pumpinfo[PumpNo,5] # q1, m^3/h
alpha2 <-Pumpinfo[PumpNo,7] # q2
alpha3 <-Pumpinfo[PumpNo,9] # q3

c <-(((beta3-beta1)-((beta2-beta1)*(alpha3-alpha1))/(alpha2-alpha1))/
      ((alpha3^2-alpha1^2)-((alpha2^2-alpha1^2)*(alpha3-alpha1))/(alpha2-alpha1)))
b <-((beta2-beta1)-c*(alpha2^2-alpha1^2))/(alpha2-alpha1)
a <-beta1-(b*alpha1+c*alpha1^2)

# test the gernealised pump curve.

q <- seq(0,360, length=25)
h <- a+b*q+c*q^2

library(ggplot2)
qplot(q,h, geom="line",main="Base Pump Curve", xlab="Nm^3/h", ylab="meters WG")

df<-data.frame(q=q, h=h)
ggplot(df, aes(x=q,y=h))+
  geom_line(data=df,lwd=0.75)+
  geom_point(data=df[1,],colour='red',size=5)+
  geom_point(data=df[19,],colour='red',size=5)+
  geom_point(data=df[25,],colour='red',size=5)


# Calculate benchmark data set for RPM vs Flow, Pressure, Power assuming a constant wheel diameter.
# where dp1 / dp2 = (RPM1 / RPM2)^2   or dp2 = dp1 / (RPM1 /RPM2)^2

RPM <- seq(from=2450,to=3050,length=25) # Range of RPM rates
PrHead=outer(RPM,q,function(RPM,q) (a+b*q+c*q^2)*(RPM/RPM1)^2)
# head(PrHead,5)


# contour (RPM,q,PrHead, xlab="x, m^3/h", ylab="y, meters WG")
# image(RPM,q,PrHead, xlab="x, m^3/h", ylab="y, meters WG")
# persp(RPM,q,PrHead, xlab="x, m^3/h", ylab="y, meters WG")
# persp(RPM,q,PrHead ,theta =30, phi =20, xlab="x, m^3/h", ylab="y, meters WG")

MinWPS <- 2000 # minimum water pump speed
MinZ <- 85


# read the actual pump performance data and convert into units similar to the theoretical pump curve.
actx1<-TOT111[,38]  # Water Pump Actual RPM
acty1<-TOT111[,42]*60/1000  # Actual flow rate, litres/min -> m^3/h
actz1<-TOT111[,24]*10  # convert barg to mWG
actdf1<-data.frame(actx1,acty1,actz1)  # actual data set converted to a dataframe
actwp1<-subset(actdf1,actx1>MinWPS & actz1>MinZ)  # remove null points

actx2<-TOT112[,38]  # Water Pump Actual RPM
acty2<-TOT112[,42]*60/1000  # Actual flow rate, litres/min -> m^3/h
actz2<-TOT112[,24]*10  # convert barg to mWG
actdf2<-data.frame(actx2,acty2,actz2)  # actual data set converted to a dataframe
actwp2<-subset(actdf2,actx2>MinWPS & actz2>MinZ)  # remove null points

actx3<-MAT102[,96]  # Water Pump Actual RPM
acty3<-MAT102[,100]*60/1000  # Actual flow rate, litres/min -> m^3/h
actz3<-MAT102[,82]*10  # convert barg to mWG
actdf3<-data.frame(actx3,acty3,actz3)  # actual data set converted to a dataframe
actwp3<-subset(actdf3,actx3>MinWPS & actz3>MinZ)  # remove null points



# Lollipop Evaluation
library(scatterplot3d)
library(rgl)

lollipop3d <- function(data.x,data.y,data.z,surf.fun,surf.n=50,
                       xlim=range(data.x),
                       ylim=range(data.y),
                       zlim=range(data.z),
                       asp=c(y=1,z=1),
                       xlab=deparse(substitute("RPM")),
                       ylab=deparse(substitute("m^3/h")),
                       zlab=deparse(substitute("mWG")),
                       alpha.surf=0.4,
                       col.surf=fg,col.stem=c(fg,fg),
                       col.pt="gray",type.surf="line",ptsize,
                       lwd.stem=2,lit=TRUE,bg="white",fg="black",
                       col.axes=fg,col.axlabs=fg,
                       axis.arrow=TRUE,axis.labels=TRUE,
                       box.col=bg,
                       axes=c("lines","box")) {
  axes <- match.arg(axes)
  col.stem <- rep(col.stem,length=2)
  x.ticks <- pretty(xlim)
  x.ticks <- x.ticks[x.ticks>=min(xlim) & x.ticks<=max(xlim)]
  x.ticklabs <- if (axis.labels) as.character(x.ticks) else NULL
  y.ticks <- pretty(ylim)
  y.ticks <- y.ticks[y.ticks>=min(ylim) & y.ticks<=max(ylim)]
  y.ticklabs <- if (axis.labels) as.character(y.ticks) else NULL
  z.ticks <- pretty(zlim)
  z.ticks <- z.ticks[z.ticks>=min(zlim) & z.ticks<=max(zlim)]
  z.ticklabs <- if (axis.labels) as.character(z.ticks) else NULL
  if (!missing(surf.fun)) {
    surf.x <- seq(xlim[1],xlim[2],length=surf.n)
    surf.y <- seq(ylim[1],ylim[2],length=surf.n)
    surf.z <- outer(surf.x,surf.y,surf.fun)  ## requires surf.fun be vectorized
    z.interc <- surf.fun(data.x,data.y)
    zdiff <- diff(range(c(surf.z,data.z)))
  } else {
    z.interc <- rep(min(data.z),length(data.x))
    zdiff <- diff(range(data.z))
  }
  xdiff <- diff(xlim)
  ydiff <- diff(ylim)
  y.adj <- if (asp[1]<=0) 1 else asp[1]*xdiff/ydiff
  data.y <- y.adj*data.y
  y.ticks <- y.adj*y.ticks
  ylim <- ylim*y.adj
  ydiff <- diff(ylim)
  z.adj <- if (asp[2]<=0) 1 else asp[2]*xdiff/zdiff
  data.z <- z.adj*data.z
  if (!missing(surf.fun)) {
    surf.y <- y.adj*surf.y
    surf.z <- z.adj*surf.z
  }
  z.interc <- z.adj*z.interc
  z.ticks <- z.adj*z.ticks
  zlim <- z.adj*zlim
  open3d()
  clear3d("all")
  light3d()
  bg3d(color=c(bg,fg))
  if (!missing(surf.fun)) 
    surface3d(surf.x,surf.y,surf.z,alpha=alpha.surf,
              front=type.surf,back=type.surf,
              col=col.surf,lit=lit)
  if (missing(ptsize)) ptsize <- 0.01*xdiff
  ## draw points
  spheres3d(data.x,data.y,data.z,r=ptsize,lit=lit,color=col.pt)
## draw lollipops
#  apply(cbind(data.x,data.y,data.z,z.interc),1,
#        function(X) {
#          lines3d(x=rep(X[1],2),
#                  y=rep(X[2],2),
#                  z=c(X[3],X[4]),
#                  col=ifelse(X[3]>X[4],col.stem[1],
#                             col.stem[2]),lwd=lwd.stem)
#        })
  bbox <- par3d("bbox")
  if (axes=="box") {
    bbox3d(xat=x.ticks,xlab=x.ticklabs,
           yat=y.ticks,ylab=y.ticklabs,
           zat=z.ticks,zlab=z.ticklabs,lit=lit)
  } else if (axes=="lines") { ## set up axis lines
    bbox <- par3d("bbox")
    axis3d(edge="x",at=x.ticks,labels=x.ticklabs,
           col=col.axes,arrow=axis.arrow)
    axis3d(edge="y",at=y.ticks,labels=y.ticklabs,
           col=col.axes,arrow=axis.arrow)
    axis3d(edge="z",at=z.ticks,labels=z.ticklabs,
           col=col.axes,arrow=axis.arrow)
    box3d(col=col.axes)
  }
  decorate3d(xlab=xlab, ylab=ylab, zlab=zlab, box=FALSE, axes=FALSE, col=col.axlabs)
}

## lollipops only ...Note need to add back the code that draws the lollipops (lines 139 to 146)


lollipop3d(actwp1$actx1,actwp1$acty1,actwp1$actz1)
lollipop3d(actwp2$actx2,actwp2$acty2,actwp2$actz2)
lollipop3d(actwp3$actx3,actwp3$acty3,actwp3$actz3)


# Quick plot in 3d showing 3 data sets in a common frame. 
plot3d(actwp1$actx1,actwp1$acty1,actwp1$actz1, r=5,col='red' )
spheres3d(actwp2$actx2,actwp2$acty2,actwp2$actz2, r=2, color="blue")
spheres3d(actwp3$actx3,actwp3$acty3,actwp3$actz3, r=2, color="green")

# setwd("C:/Users/warwi/OneDrive/R Language/waterpump3D")
filename <- writeWebGL(dir = file.path(getwd(), "webGL"), 
                       width = 500, reuse = TRUE)
attr(filename, "reuse")
if (interactive())
  browseURL(paste0("file://", filename))


## lollipops plus theoretical surface
PrFun <- function(RPM,q) (RPM/RPM1)^2*(a+b*q+c*q^2)  
lollipop3d(actwp1$actx1,actwp1$acty1,actwp1$actz1,PrFun,col.pt="red",col.stem=c("red","blue"))



