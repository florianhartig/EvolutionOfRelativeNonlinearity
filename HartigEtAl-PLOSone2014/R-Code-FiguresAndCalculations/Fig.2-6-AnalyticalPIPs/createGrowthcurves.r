################################################################################
#
# script for visualizing density dependence - results are not used in the article
#
# http://florianhartig.wordpress.com/
################################################################################


rm(list=ls(all=TRUE)) 
setwd("/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.2-AnalyticalPIPs/")
source("functions.R")



pdf("growthcurves.pdf", width = 8, height = 8)



par(mfrow=c(2,2))
color = gray(0:30/30)
x=seq(0,2,by=0.02)
bval <- 10^(seq(-0.5,1,length.out=30))

plot(x, mss(x, b=6), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="Original MSS",  lty = 3 , lwd = 2, col = color[1])
for (i in 1:length(bval)){
  lines(x, mss(x, b=bval[i]) , col = color[i],  lty = 1, lwd = 1.5)  
}

lines(x, mss(x, b=bval[9]) , col = "darkred",  lty = 2, lwd = 2)  
lines(x, mss(x, b=bval[24]) , col = "darkgreen",  lty = 2, lwd = 2)  

abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)

plot(x, mssmod(x, b=6), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="Modified MSS",  lty = 3 , lwd = 2, col = color[1])
for (i in 1:length(bval)){
  lines(x, mssmod(x, b=bval[i]) , col = color[i],  lty = 1, lwd = 1.5)  
}

lines(x, mssmod(x, b=bval[9]) , col = "darkred",  lty = 2, lwd = 2)  
lines(x, mssmod(x, b=bval[24]) , col = "darkgreen",  lty = 2, lwd = 2)  

abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)







x=seq(0.8,1.2,by=0.02)


plot(x, log(mss(x, b=6)), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="",  lty = 3 , lwd = 2, col = color[1])
for (i in 1:length(bval)){
  lines(x, log(mss(x, b=bval[i])) , col = color[i],  lty = 1, lwd = 1.5)  
}
abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)



plot(x, log(mssmod(x, b=6)), type = "n",xlab = "N/K", ylab = "Reproductive rate", main="",  lty = 3 , lwd = 2, col = color[1])
for (i in 1:length(bval)){
  lines(x, log(mssmod(x, b=bval[i])) , col = color[i],  lty = 1, lwd = 1.5)  
}
abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
legend("topright", bg="white", cex = 0.7,inset=0.0, legend = round(bval, digits = 1), lwd = 1.5, col = color, merge = TRUE)


dev.off()



