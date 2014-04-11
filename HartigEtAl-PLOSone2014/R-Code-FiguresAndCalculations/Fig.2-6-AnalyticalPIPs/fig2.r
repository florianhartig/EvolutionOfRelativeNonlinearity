################################################################################
#
#
# http://florianhartig.wordpress.com/
################################################################################


rm(list=ls(all=TRUE)) 
setwd("/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/R-Code-FiguresAndCalculations/Fig.2-6-AnalyticalPIPs")
source("functions.R")

require(fields)
require(graphics)


#testresponse<- function(n)mss(n, r=5, b=1, K=1, d=0.0, procerror = T, max = 0.03)
#plot_popdyn(create_resident_dynamics(testresponse))

invasibility_MSS = create_PIP(popfun = mss)
invasibility_Ricker = create_PIP(popfun = ricker)
invasibility_Hassell = create_PIP(popfun = hassell)


globalmin = min(c(min(invasibility_MSS[[2]]),min(invasibility_Ricker[[2]]), min(invasibility_Hassell[[2]]) )) 
globalmax = max(c(max(invasibility_MSS[[2]]),max(invasibility_Ricker[[2]]), max(invasibility_Hassell[[2]]) )) 
zlimits = c(globalmin, globalmax)
        
col = rescale_colors(basecol, 0.75, 3, 3)



################################################################################
#
#  Plots for the inlets


png("temp1.png")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))

image(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),ifelse((invasibility_MSS[[2]]>0) * (t(invasibility_MSS[[2]]) > 0) ,1,NA) , col = "gray22", axes = F, 
      main = "", xlab = "",ylab = "")

bcri = log10(2.51)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
dev.off()



png("temp2.png")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))

image(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),ifelse((invasibility_Ricker[[2]]>0) * (t(invasibility_Ricker[[2]]) > 0) ,1,NA) , col = "gray22", axes = F, 
      main = "", xlab = "",ylab = "")

bcri = log10(4.45725)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
dev.off()



png("temp3.png")

par(mar=c(0,0,0,0), oma=c(0,0,0,0))

image(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),ifelse((invasibility_Hassell[[2]]>0) * (t(invasibility_Hassell[[2]]) > 0) ,1,NA) , col = "gray22", axes = F, 
      main = "", xlab = "",ylab = "")

bcri = log10(2.64223)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
dev.off()





################################################################################
#
#  Main plots


qual = 5

#tiff("fig2.tiff", width= 1500 * qual, height = 520 * qual, res = 72 * qual, compression = "zip")
pdf("fig2-raw.pdf", width = 20, height = 7)
#graphics.off()
#windows(width= 1500, height = 520)

scaling = 2
par(mfrow= c(1,3))
par(mar = c(3*scaling,3*scaling,3*scaling,3*scaling))





labcex = scaling
maincex = scaling
par(cex.axis = scaling)
plussize = 1.5*scaling
mpg = c(2*scaling, 1,0)



image(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),invasibility_MSS[[2]], axes = F, zlim  = zlimits,
      main = "(a) MSS", font.main = 1, xlab = "Resident density-compensation strategy b",ylab = "Invading density-compensation strategy b",
      col = col , cex.lab = labcex, cex.main= maincex, mgp = mpg)
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),invasibility_MSS[[2]], legend.only = T, zlim  = zlimits, col = col)
#contour(log10(bvalues), log10(bvalues), invasibility_dist_cat, add = T, nlevels = 1, lty = 3, lwd = 2, drawlabels = F)
image(log10(invasibility_MSS[[1]]),log10(invasibility_MSS[[1]]),ifelse(invasibility_MSS[[2]]<0.00001,1,NA), add =T, col = "#00000030")
bcri = log10(2.51)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
#abline(h = bcri, lwd = 1, lty = 2)
text(log10(1.1), log10(0.4), "-", cex = plussize)
text(log10(5), log10(11), "-", cex = plussize)
text(log10(0.8), bcri, "+", cex = plussize)
text(log10(8), bcri, "+", cex = plussize)


image(log10(invasibility_Ricker[[1]]),log10(invasibility_Ricker[[1]]),invasibility_Ricker[[2]], axes = F, zlim  = zlimits,
      main = "(b) Ricker", font.main = 1, xlab = "Resident density-compensation strategy b", ylab = "", col = col, cex.lab = labcex, cex.main= maincex, mgp = mpg)
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(log10(invasibility_Ricker[[1]]),log10(invasibility_Ricker[[1]]),invasibility_Ricker[[2]], legend.only = T, zlim  = zlimits, col = col)
#contour(log10(bvalues_ricker), log10(bvalues_ricker), invasibility_ricker_cat, add = T, nlevels = 1, lty = 3, lwd = 2, drawlabels = F)
image(log10(invasibility_Ricker[[1]]),log10(invasibility_Ricker[[1]]),ifelse(invasibility_Ricker[[2]]<0.00001,1,NA), add =T, col = "#00000030")
bcri = log10(4.45725)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
#abline(h = bcri, lwd = 1, lty = 2)
text(log10(2), log10(0.4), "-", cex = plussize)
text(log10(6), log10(11), "-", cex = plussize)
text(log10(1.8), bcri, "+", cex = plussize)
text(log10(10), bcri, "+", cex = plussize)



image(log10(invasibility_Hassell[[1]]),log10(invasibility_Hassell[[1]]),invasibility_Hassell[[2]], axes = F, zlim  = zlimits,
      main = "(c) Hassel", font.main = 1, xlab = "Resident density-compensation strategy b", ylab = "", col = col, cex.lab = labcex, cex.main= maincex, mgp = mpg)
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(log10(invasibility_Hassell[[1]]),log10(invasibility_Hassell[[1]]),invasibility_Hassell[[2]], legend.only = T, zlim  = zlimits, col = col)
#contour(log10(bvalues_ricker), log10(bvalues_ricker), invasibility_hassel_cat, add = T, nlevels = 1, lty = 3, lwd = 2, drawlabels = F)
image(log10(invasibility_Hassell[[1]]),log10(invasibility_Hassell[[1]]),ifelse(invasibility_Hassell[[2]]<0.00001,1,NA), add =T, col = "#00000030")
#contour(log10(bvalues_ricker), log10(bvalues_ricker), invasibility_hassel_cat, add = T, nlevels = 1, lty = 3, lwd = 2, drawlabels = F)

bcri = log10(2.64223)
abline(v = bcri, lwd = 1, lty = 2)
lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
#abline(h = bcri, lwd = 1, lty = 2)
text(log10(1.1), log10(0.4), "-", cex = plussize)
text(log10(5), log10(11), "-", cex = plussize)
text(log10(0.8), bcri, "+", cex = plussize)
text(log10(8), bcri, "+", cex = plussize)


dev.off()


