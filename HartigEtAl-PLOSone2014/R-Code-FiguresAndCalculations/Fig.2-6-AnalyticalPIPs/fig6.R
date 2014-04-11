################################################################################
#
# http://florianhartig.wordpress.com/
################################################################################


rm(list=ls(all=TRUE)) 
setwd("/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.2-AnalyticalPIPs")
source("functions.R")

require(fields)
require(graphics)

invasibility_MSSmod = create_PIP(popfun = mssmod)




graphics.off()
#png(width= 1000, height = 340, res = 300)

pdf("fig6-raw.pdf", width = 20, height = 7)

scaling = 2

#png("fig6.png", width= 1500, height = 520)

par(mfrow=c(1,3), mar = c(3*scaling,3*scaling,3*scaling,3*scaling) )

labcex = scaling
maincex = scaling
par(cex.axis = scaling)
plussize = 1.5*scaling
mpg = c(2*scaling, 1,0)



color = gray(0:30/30)
x=seq(0,2,by=0.02)
bval <- 10^(seq(-0.5,1,length.out=30))

plot(x, mss(x, b=6), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="Original MSS",  lty = 3 , lwd = 2, col = color[1]
     , cex.lab = labcex, cex.main= maincex, mgp = mpg, font.main = 1)
for (i in 1:length(bval)){
  lines(x, mss(x, b=bval[i]) , col = color[i],  lty = 1, lwd = 1.5)  
}

lines(x, mss(x, b=bval[9]) , col = "darkred",  lty = 2, lwd = 2)  
lines(x, mss(x, b=bval[24]) , col = "darkgreen",  lty = 2, lwd = 2)  

abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
#legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)

plot(x, mssmod(x, b=6), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="Modified MSS",  lty = 3 , lwd = 2, col = color[1]
     , cex.lab = labcex, cex.main= maincex, mgp = mpg, font.main = 1)
for (i in 1:length(bval)){
  lines(x, mssmod(x, b=bval[i]) , col = color[i],  lty = 1, lwd = 1.5)  
}

lines(x, mssmod(x, b=bval[9]) , col = "darkred",  lty = 2, lwd = 2)  
lines(x, mssmod(x, b=bval[24]) , col = "darkgreen",  lty = 2, lwd = 2)  

abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
#legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)




globalmin = min(invasibility_MSSmod[[2]]) 
globalmax = max(invasibility_MSSmod[[2]]) 
zlimits = c(globalmin, globalmax)

col = rescale_colors(basecol, 0.67, 3, 3)


image(log10(invasibility_MSSmod[[1]]),log10(invasibility_MSSmod[[1]]),invasibility_MSSmod[[2]], axes = F, zlim  = zlimits,
      main = "Invasibility modified MSS", font.main = 1, xlab = "Resident density-compensation strategy b",ylab = "Invading density-compensation strategy b",
      col = col , cex.lab = labcex, cex.main= maincex, mgp = mpg)
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(log10(invasibility_MSSmod[[1]]),log10(invasibility_MSSmod[[1]]),invasibility_MSSmod[[2]], legend.only = T, zlim  = zlimits, col = col)
#contour(log10(invasibility_MSSmod[[1]]),log10(invasibility_MSSmod[[1]]),invasibility_MSSmod[[2]], add = T, nlevels = 1, lty = 3, lwd = 2, drawlabels = F)
image(log10(invasibility_MSSmod[[1]]),log10(invasibility_MSSmod[[1]]),ifelse(invasibility_MSSmod[[2]]<0.000001,1,NA), add =T, col = "#00000030")
bcri = log10(2.51)
#abline(v = bcri, lwd = 1, lty = 2)
#lines(c(bcri + 0.05, bcri-0.05), c(bcri,bcri))
#abline(h = bcri, lwd = 1, lty = 2)
text(log10(1.1), log10(0.4), "-", cex = plussize)
text(log10(6), log10(11), "-", cex = plussize)
text(log10(0.5), log10(4), "+", cex = plussize)
text(log10(10), log10(1), "+", cex = plussize)

#contour(log10(bvalues), log10(bvalues), mutual, add = T, levels = c(1000), lty = 2, lwd = 1, drawlabels = F)

dev.off()
