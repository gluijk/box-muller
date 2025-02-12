# Box–Muller transform to generate normally distributed random numbers
# from uniformly distributed random numbers
# www.overfitting.net
# https://www.overfitting.net/2025/02/la-transformada-de-box-muller-con-r.html


library(tiff)  # save 16-bit TIFF's
library(data.table)


rnorm.boxmuller=function(U, mean=0, sd=1) {
    # Box–Muller transform to get normally distributed pairs
    # from a M x 2 matrix of [0,1] uniformly distributed pairs
    R=sqrt(-2*log(U[,1]))
    theta=2*pi*U[,2]
    Z0=R*cos(theta)
    Z1=R*sin(theta)
    return(cbind(Z0,Z1)*sd + mean)
}

NewBitmap = function(dimx, dimy, val=0) {
    # Crea bitmap de dimensiones dimx y dimy
    return(array(val,c(dimx,dimy)))
}

# Por Carlos Gil Bellosta
indices.drawline = function(x0, y0, x1, y1) {
    x0=round(x0)
    x1=round(x1)
    y0=round(y0)
    y1=round(y1)
    
    if (y0 == y1) return(cbind(x0:x1, y0)) # Recta de m=0 o un punto
    if (abs(x1 - x0) >= abs(y1 - y0)) { # Recta de 0 < |m| <= 1
        m = (y1 - y0) / (x1 - x0)
        cbind(x0:x1, round(y0 + m * ((x0:x1) - x0)))
    } else indices.drawline(y0, x0, y1, x1)[, 2:1]  # Recta de |m| > 1
    # Llamada traspuesta recursiva y traspuesta
}

DrawLine = function(img, x0, y0, x1, y1, inc=TRUE, val=1) {
    # Dibuja recta desde (x0,y0)-(x1,y1)
    # Por defecto método no destructivo y con valor=1
    indices=indices.drawline(x0, y0, x1, y1)
    if (inc) img[indices]=img[indices]+val
    else img[indices]=val
    
    return(img)
}

SaveBitmap = function(img, name, trunc=TRUE, gamma=1) {
    # Guarda bitmap en formato PNG
    # Solo si trunc=FALSE y la imagen excede de 1 se reescala a 1
    require(tiff)
    img[img<0]=0
    if (trunc) img[img>1]=1
    if (tolower(substr(name, nchar(name)-3, nchar(name))) != ".tif") name=paste0(name,".tif")
    # writePNG(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name)
    writeTIFF(t(img[,ncol(img):1] / max(max(img),1))^(1/gamma), name,
              bits.per.sample=16)
}


##################################

# Scatterplots
N=6000
datosunif=matrix(runif(N), ncol=2)
datosnorm=rnorm.boxmuller(datosunif)

dev.off()
par(mfrow = c(2, 1))
par(cex=1)
plot(datosunif, col=rgb(1, 0, 0, alpha=0.04), pch=16, cex=1, asp=1,
     xlim=c(-0.5,1.5), ylim=c(-0.5,1.5), xlab='U1', ylab='U2',
     xaxt="n", yaxt="n")
axis(side=1, at=seq(-1,2,1))
axis(side=2, at=seq(-1,2,1))
abline(h=0, v=0)
U1=0.1  # example
U2=0.8
points(U1, U2, pch=16, cex=1.5)

plot(datosnorm, col=rgb(1, 0, 0, alpha=0.04), pch=16, cex=1, asp=1,
     xlim=c(-3,3), ylim=c(-3,3))
axis(side=1, at=seq(-3,3,1))
axis(side=2, at=seq(-3,3,1))
abline(h=0, v=0)
symbols(0, 0, circles=1, add=TRUE, inches=FALSE)
example=rnorm.boxmuller(matrix(c(U1,U2),ncol=2))
points(example[1], example[2], pch=16, cex=1.5)


# Histograms
N=4000000
datosunif=matrix(runif(N), ncol=2)
datosnorm=rnorm.boxmuller(datosunif)

dev.off()
par(mfrow = c(2, 1))
hist(datosunif, xlim=c(-0.5,1.5), breaks=800)
hist(datosnorm, xlim=c(-4,4), breaks=800)


# Transformation
NPLOT=10000
dev.off()
plot(0, 0, xlim=c(-3, 3), ylim=c(-3, 3), type="n", xlab="X", ylab="Y", asp=1)
abline(h=c(0,1), v=c(0,1), col='lightgray')
symbols(0, 0, circles=1, fg='lightgray', add=TRUE, inches=FALSE)
axis(side=1, at=seq(-3,3,1))
axis(side=2, at=seq(-3,3,1))
segments(datosunif[(1:NPLOT), 1], datosunif[(1:NPLOT), 2],
         datosnorm[(1:NPLOT), 1], datosnorm[(1:NPLOT), 2],
         col=rgb(1, 0, 0, alpha=0.01), lwd=1)


# Hires transformation
N=800000
datosunif=matrix(runif(N), ncol=2)
datosnorm=rnorm.boxmuller(datosunif)

DT=as.data.table(cbind(datosunif, datosnorm))
colnames(DT)=c("U1", "U2", "Z0", "Z1")  # Wikipedia naming

DIMY=1080
DIMX=1080
OFFY=DIMY/2
OFFX=DIMX/2

img=NewBitmap(DIMY, DIMY)
MAXX=3  # must be >=1
MAXY=3  # must be >=1

DT$U1=round(DT$U1*OFFX/MAXX + OFFX)
DT$U2=round(DT$U2*OFFY/MAXY + OFFY)
DT$Z0=round(DT$Z0*OFFX/MAXX + OFFX)
DT$Z1=round(DT$Z1*OFFY/MAXY + OFFY)

# Drop out of range values and group repeated values
DT=DT[Z0 >= 1 & Z0 <= DIMX & Z1 >=1 & Z1 <= DIMY, .(count = .N),
       by=.(U1, U2, Z0, Z1)]
M=as.matrix(DT)
rm(datosunif)
rm(datosnorm)
rm(DT)

U1=M[,1]
U2=M[,2]
Z0=M[,3]
Z1=M[,4]
val=M[,5]

for (i in 1:nrow(M)) img=DrawLine(img, U1[i], U2[i], Z0[i], Z1[i], val[i])
SaveBitmap(img, "boxmuller.tif", trunc=FALSE, gamma=2.2)

# Generate separate axes
img=NewBitmap(DIMY, DIMY)
img=DrawLine(img, 1, OFFY, DIMX, OFFY)
img=DrawLine(img, OFFX, 1, OFFX, DIMY)
SaveBitmap(img, "boxmulleraxes.tif")

