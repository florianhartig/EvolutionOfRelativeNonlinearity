################################################################################
#
#
# http://florianhartig.wordpress.com/
################################################################################

rm(list=ls(all=TRUE)) 

#install.packages("tgp")
require(fields)
require(tgp)
require(colorRamps)



wdir <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.3-IBM_InvasionAndCoexistence"
setwd(wdir)

survival <- read.csv("data/Coexistence competitive-exclusion-table.csv", skip = 6, )

survival = survival[order(survival$X.run.number.),] 
#summary(survival)

data = survival

x1 <- unique(data$initial.strategy.1)
x2 <- unique(data$initial.strategy.2)

species1 <- (matrix(as.numeric(data$X.count.species1), byrow=T, ncol=length(x1))) - 1 
species2 <- ((matrix(as.numeric(data$X.count.species2), byrow=T, ncol=length(x1))) - 1 )  * 2 
ticks <- log10(matrix(as.numeric(data$ticks), byrow=T, ncol=length(x1)) )

species12 <- species1 + species2  
species12[species12 == 0] = NA

species12vector <- as.numeric(data$X.count.species1) + 2* as.numeric(data$X.count.species2) -3
species12vector[species12vector == 0] = NA

numberofruns = 36    # number of runs for which invasion was tested
invasion_all <- read.csv("data/Coexistence invasion-stochastic-all-table.csv", skip = 6, )
invasion_all = invasion_all[order(invasion_all$X.run.number.),]
invasibility_all <- (matrix(as.numeric(invasion_all$invasioncounter ) , byrow=T, ncol=length(x1))) 




#invasion_smoothened = interp.loess(log10(invasion_all$initial.strategy.1), log10(invasion_all$initial.strategy.2), 1/invasion_all$invasioncounter, gridlen = c(251,251))
dominance_smoothened = interp.loess(log10(data$initial.strategy.1), log10(data$initial.strategy.2), species12vector, gridlen = c(251,251), span = 0.02)


mutual = interp.loess(log10(data$initial.strategy.1),log10(data$initial.strategy.2), as.vector((1-invasibility_all/numberofruns) * t(1-invasibility_all/numberofruns)), gridlen = c(251,251),  span = 0.01)
mutual = as.matrix(mutual$z)
mutual[mutual < 0.1] = NA



#invasion_all_indep = cbind(invasion_all$initial.strategy.1,invasion_all$initial.strategy.2)
#fit<- Tps(invasion_all_indep, invasion_all$invasioncounter)
           


colorscheme = colorRampPalette(c("darkgreen", "yellow", "grey99")) 
colorscheme2 = colorRampPalette(c("white","lightgrey", "black"))
colorscheme3 = colorRampPalette(c("lightgrey", "black"))
colorscheme4 = colorRampPalette(c("white", "grey90", "grey75", "grey60", "grey45", "grey30", "grey15", "grey9", "black"))


#graphics.off()
#windows(width= 600, height = 310)
pdf("fig3-raw.pdf", width = 20, height = 10)
oldpar <- par(mfrow=c(1,2)) # set graphic device

#image.plot(log(x1), log(x2), species12, col = colorscheme2(3), main = "Species dominance", xlab = "ln b species 1 ", ylab = "ln b species 2",  nlevel = 3, legend.mar = 2)      #



#par(oma = c(1,1,3,1))
#mtext("Evolutionary vs. dynamic stability\n", outer=T, cex = 1.8)

#1-(1-1/invasibility_all)^(1/3)

par(mar = c(4.5,4.5,4.5,4.5), cex = 1.4, cex.lab = 1.1)
image(log10(x1), log10(x2), 1 - (invasibility_all/numberofruns)^(1/3),  col = colorscheme4(100),  axes = F, main = "(a) Invasion probability", font.main = 1, xlab = " Density-compensation strategy, b, of resident", ylab = " Density-compensation strategy, b, of mutant")
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(1 - (invasibility_all/numberofruns)^(1/3), legend.only = T,  col = colorscheme4(100))
abline(v = log10(2.51), lwd = 2, col = "darkgrey" , lty = 2)
abline(h = log10(2.51), lwd = 2, col = "darkgrey",  lty = 2)

image(log10(x1), log10(x2), mutual, add =T, col = "#FF000050")

#contour(mutual, add = T, levels =c(0.1) , lwd = 3, lty = 3,  drawlabels = F, col = "red")

#text(log10(0.9), log10(6), "MI", cex = 1.6, col = "red")
#text(log10(6), log10(0.9), "MI", cex = 1.6, col = "red")


par(mar = c(4.5,4.5,4.5,4.5), cex = 1.4, cex.lab = 1.1)

zlimits = c(0,6.7)
labelvalues = 10^seq(0,6.7, length = 65)
image(log10(x1), log10(x2), ticks , col = tim.colors(64), main = "(b) Time until competitive exclusion", font.main = 1, axes = F, xlab = " Density-compensation strategy, b, of species 1", ylab = " Density-compensation strategy, b, of species 2")
axis(side = 1, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
axis(side = 2, at = log10(c(seq(from = 0.1, to = 1, by = 0.1), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F )
axis(side = 2, at = log10(c(0.3, 1, 3, 10)), labels = c("0.3", "1", "3", "10"), lwd.ticks = 2, tcl = - 0.5 )
image.plot(ticks, legend.only = T, lab.breaks = as.character(round(labelvalues)))
contour(dominance_smoothened,  add = T, levels =c(1.5) , lwd = 2, lty = 2,  drawlabels = F)
text(log10(0.4), log10(2), "2", cex = 2)
text(log10(10), log10(2), "2", cex = 2)
text(log10(2.4), log10(0.5), "1", cex = 2)
text(log10(2.4), log10(10), "1", cex = 2)
contour(mutual, add = T, levels =c(0.1) , lwd = 3, lty = 3,  drawlabels = F, col = "red")



#
#
#par(mar = c(4.5,4.5,4.5,4.5), cex = 1.2)
#image.plot(log10(x1), log10(x2), 1/invasibility_all,  col = colorscheme4(100), main = "(a) Invasion probability", xlab = "Log b resident ", ylab = "Log b invader")  
#abline(v = log10(2.51), lwd = 2, col = "darkgrey" , lty = 2)
#abline(h = log10(2.51), lwd = 2, col = "darkgrey",  lty = 2)
##points(log(2.51),log(2.51), col="red", cex = 5, lwd = 4  )
#
#
#
#par(mar = c(6,6,6,6),cex = 1.2)
#image.plot(log10(x1), log10(x2), ticks ,  main = "(b) Log time to\ncopetitive exclusion", xlab = "Log b species I ", ylab = "Log b species II", zlim = c(0,6.7), lab.breaks = c("1", "10", "10^2", "10^3", "10^4", "10^5", "10^6"))
#
##contour(invasion_smoothened, add = T, levels =c(0.1) , lwd = 2, drawlabels = F)
#contour(dominance_smoothened, add = T, levels =c(1.5) , lwd = 2, lty = 2,  drawlabels = F)
#text(-0.5, 0.5, "sp. 2", cex = 2)
#text(1.5, 0.8, "sp. 2", cex = 2)
#text(0.8, -0.8, "sp. 1", cex = 2)
#text(0.8, 1.5, "sp. 2", cex = 2)
#
#

#legend("bottomright", bg="white", inset=0.0, legend = c("species 1", "species2"), col = c("red", "blue"), merge = TRUE)

dev.off()
par <- oldpar

