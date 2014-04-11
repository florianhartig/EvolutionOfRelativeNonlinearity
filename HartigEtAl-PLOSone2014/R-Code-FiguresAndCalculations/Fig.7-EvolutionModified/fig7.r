################################################################################
#
# http://florianhartig.wordpress.com
################################################################################

rm(list=ls(all=TRUE)) 

library(fields)
require(tgp)
require(colorRamps)

################################################################################
# read in all data 

setwd("/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.7-EvolutionModified")

timesteps = 100
totaltime = 150000
lines = totaltime / timesteps + 1

bval = seq(0.025,6.975, length = 140)
timeval = seq(0,150000, by=100)  / 1000

# square brackets need to be removed before reading in the data

evoldata = read.table("rate-5_disturbance-0.05startingpoint-4.dat", header = T, na.strings = "0")

evoldata = as.matrix(evoldata[1:lines,])

################################################################################


scaling = 1.7


#graphics.off()
#windows(width= 1000, height = 520)

pdf("fig7-raw.pdf", width = 20, height = 10)

labcex = scaling
maincex = scaling
par(cex.axis = scaling)
plussize = 1.5*scaling
mpg = c(2*scaling, 1,0)

par(mfrow = c(1,2), mar = c(3*scaling,3*scaling,3*scaling,3*scaling))


par(cex.lab = 1.2)
pal = colorRampPalette(c("orange", "yellow", "white"), bias = 1)
image(timeval, bval, evoldata , col = pal(100), main = "(A) Evolutionary branching", font.main=1, xlab = "Time (10^3 generations)", ylab = "Density compensation strategy b"
      , cex.lab = labcex, cex.main= maincex, mgp = mpg, font.main = 1)
abline(v = 30)



########################################################################################

survival <- read.csv("Coexistence competitive-exclusion-MSS-Mod-table.csv", skip = 6, )

survival = survival[order(survival$X.run.number.),] 

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

dominance_smoothened = interp.loess(log10(data$initial.strategy.1), log10(data$initial.strategy.2), species12vector, gridlen = c(251,251), span = 0.02)


zlimits = c(0,6.7)
labelvalues = 10^seq(0,6.7, length = 65)
image(log10(x1), log10(x2), ticks , col = tim.colors(64), main = "(b) Time until competitive exclusion", axes = F, xlab = " Density-compensation strategy, b, of species 1", ylab = " Density-compensation strategy, b, of species 2"
      , cex.lab = labcex, cex.main= maincex, mgp = mpg, font.main = 1)
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


dev.off()
