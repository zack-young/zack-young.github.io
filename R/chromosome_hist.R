#/bin/Rscript

plotChrom <- function(xleft, ybottom, height, len, centro, binsize){
  
  len <- len/binsize
  centro <- centro/binsize
  
  # r vertical, r horizontal
  rv <- height/2
  rh <- len/70
  rs <- seq(0,pi,len=100)
  
  # left semi-circle
  lx <- c(xleft, xleft + rh - rh*sin(rs), xleft, xleft)
  ly <- c(ybottom, ybottom + rv - rv*cos(rs), ybottom + height, ybottom)
  polygon(lx, ly, border = "white", col = "white")
  
  # right semi-circle
  rx <- c(xleft + len,xleft + len, xleft + len - rh + rh*sin(rs), xleft + len)
  ry <- c(ybottom, ybottom + height, ybottom + rv +rv*cos(rs), ybottom)
  polygon(rx, ry, border = "white", col = "white")
  
  # top tri
  tx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  ty <- c(ybottom + height, ybottom + rv, ybottom + height)
  polygon(tx, ty, border = "white", col = "white")
  
  # bottom tri
  bx <- c(centro - 0.5*rh, centro, centro + 0.5*rh)
  by <- c(ybottom, ybottom + rv, ybottom)
  polygon(bx, by, border = "white", col = "white")
  
  # chrom border
  if (is.na(centro)){
    cbx <- c(xleft + rh - rh*sin(rs), xleft + len - rh + rh*sin(rs))
    cby <- c(ybottom + rv - rv*cos(rs), ybottom + rv + rv*cos(rs))
    polygon(cbx, cby, border = "grey", lwd = 3)
  }else{
    cbx <- c(xleft + rh - rh*sin(rs), centro - 0.5*rh, centro, centro + 0.5*rh, xleft + len - rh + rh*sin(rs), centro + 0.5*rh, centro, centro - 0.5*rh)
    cby <- c(ybottom + rv - rv*cos(rs), ybottom + height, ybottom + rv, ybottom + height, ybottom + rv + rv*cos(rs), ybottom, ybottom + rv, ybottom)
    polygon(cbx, cby, border = "grey", lwd = 3)
  }
}

df1 <- read.table("CS_abs_gene_sum.txt", header = F, stringsAsFactors = F)
df2 <- read.table("Zang1817_abs_gene_sum.txt", header = F, stringsAsFactors = F)

df1$V1 <- df1$V1/max(df1$V1)
df2$V1 <- df2$V1/max(df2$V1)

sarm <- 20
larm <- 30
tarm <- sarm + larm

pdf("RD_quantile.pdf", width = 10, height = 6)
par(mgp = c(1,0,0))
plot(x=NULL,y=NULL,xlim=c(-sarm,larm),ylim=c(-1.3,1.3), frame.plot = F, yaxt = "n", xaxt = "n", xlab = "Relative chromosome position", ylab = "Relative ratio", cex.lab=1.5, main = "Relative gene absence ratio between CS and Zang1817", cex.main = 1.5)
plotChrom(xleft = -sarm,ybottom = -0.05,height = 0.1,len = tarm,centro = 0.5,binsize = 1)
rect(xleft = c(-sarm:(larm-1)), ybottom = rep(0,tarm)+0.08, xright = c((-sarm+1):larm), ytop = df1$V1+0.08, col = rgb(51, 126, 187, maxColorValue=255), border = F)
# rect(xleft = seq(0,13.5,1.5), ybottom = rep(0,10)+0.03, xright = seq(1.5,15,1.5), ytop = df1$V1[11:20]+0.03, col = rgb(51, 126, 187, maxColorValue=255), border = F)
asm <- smooth.spline( seq(-sarm+0.5,larm-0.5,1), df1$V1+0.03, spar=0.3, all.knots = F, nknots = 25)
asm$y <- asm$y + 0.05
lines(asm, col="#1D4669", lwd=2, lty=2)
rect(xleft = c(-sarm:(larm-1)), ybottom = rep(0,tarm)-0.08, xright = c((-sarm+1):larm), ytop = -df2$V1-0.08, col = rgb(228, 30, 38, maxColorValue=255), border = F)
# rect(xleft = seq(0,13.5,1.5), ybottom = rep(0,10)-0.03, xright = seq(1.5,15,1.5), ytop = -df2$V1[11:20]-0.03, col = rgb(228, 30, 38, maxColorValue=255), border = F)
bsm <- smooth.spline( seq(-sarm+0.5,larm-0.5,1), -df2$V1-0.03, spar=0.3, all.knots = F, nknots = 25)
bsm$y <- bsm$y - 0.05
lines(bsm, col="#911318", lwd=2, lty=2)

text(x = -6,
     y = -0.2,
     labels = c("Absent gene in CS"),
     cex = 1,
     xpd = T,
     col = "#911318"
)

text(x = -6,
     y = 0.25,
     labels = "Absent gene in Zang1817",
     cex = 1,
     xpd = T,
     col = "#1D4669"
)


axiy = 0.08
# axis
segments(x0 = -sarm-0.5, 
         y0 = c(-axiy, axiy),
         x1 = -sarm-0.5, 
         y1 = c(-max(df2$V1)-axiy, max(df1$V1)+axiy),
         col = "black",
         xpd = T
)

ticky <- c(-max(df2$V1)-axiy, -max(df2$V1)/2 -axiy, -axiy , axiy, max(df1$V1)/2 + axiy, max(df1$V1) + axiy)

# axis ticks
segments(x0 = -sarm-0.5,
         x1 = -sarm-0.5+0.3,
         y0 = ticky,
         y1 = ticky,
         xpd = T
)

# axis texts
text(x = -sarm-0.1,
     y = ticky,
     labels = c(1, 0.5, 0, 0,0.5, 1),
     cex = 1,
     xpd = T,
     pos = 2
)

# x-axis texts
text(x = c(-14.5, 30),
     y = 0,
     labels = c("short arm","long arm"),
     cex = 1,
     xpd = T,
     pos = 2
)

dev.off()