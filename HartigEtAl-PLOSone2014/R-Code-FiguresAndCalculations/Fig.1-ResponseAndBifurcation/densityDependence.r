################################################################################
# Plot of the density response for Fig. 1 
# http://florianhartig.wordpress.com/
################################################################################

rm(list=ls(all=TRUE))

#install.packages("tgp")
require(fields)

      
wdir <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.1-ResponseAndBifurcation/"
setwd(wdir)           

mms <- function(N, r=5, b=2, K=1){
	return(r * N / (1 + (r-1) * (N/K)^b))
}

mms <- function(N, r=5, b=2, K=1, d=0.05){
	return(r * (1-d) / (1 + (r-1) * (N/K)^b) )
}



x=seq(0,2,by=0.01)
par(cex=1.2, lwd = 1.5)

color = c("#FFD700", "#FF7256", "#8B4513", "#CD5B45", "#8B3E2F" )
color = c("black", "black", "black", "gray", "gray" )
#color = rep("black",5)

plot(x, mms(x, b=6), type = "l", xlab = "", ylab = "", main="",  lty = 1 , lwd = 2, col = color[1])
#lines(x, mms(x, b=2.8), col = color[2],  lty = 1)
lines(x, mms(x, b=2.5) , col = "black",  lty = 1, lwd = 1)
lines(x, mms(x, b=1) , col = "gray",  lty = 1, lwd = 2,)
lines(x, mms(x, b=0.5) , col = "gray" ,  lty = 1, lwd = 1)
abline(v=1, lwd = 0.5, , lty = 1)
abline(h=1, lwd = 0.5, , lty = 1)
legend("topright", bg="white", col = c("black", "black", "gray", "gray"), inset=0.0, legend = c("b=6",  "b=2.5", "b=1", "b=0.1"), lwd = c(2,1,2,1), lty=1, merge = TRUE)







