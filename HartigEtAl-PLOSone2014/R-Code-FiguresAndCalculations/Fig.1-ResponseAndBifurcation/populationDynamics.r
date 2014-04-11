################################################################################
#
# florian.hartig@ufz.de       http://www.ufz.de/index.php?de=10623
################################################################################

rm(list=ls(all=TRUE))

#install.packages("tgp")
require(fields)
require(tgp)

wdir <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.1-ResponseAndBifurcation/"
setwd(wdir)

chaos <- read.csv(paste(wdir,"simulationOutput/chaos.csv", sep=""), skip = 6, )
konvergent <- read.csv(paste(wdir,"simulationOutput/konvergent.csv", sep=""), skip = 6, )
oszillierend <- read.csv(paste(wdir,"simulationOutput/oszillierend.csv", sep=""), skip = 6, )
gedaempft <- read.csv(paste(wdir,"simulationOutput/gedaempft.csv", sep=""), skip = 6, )



ltyg = "b"

#graphics.off()
#windows(width= 400, height = 400)
oldpar <- par(mfrow=c(2,2), cex.axis= 1.6, cex.lab = 1.5) # set graphic device

plot(konvergent$sum..individuals..of.turtles[1:50] / 200, type = ltyg, ylim = c(0,2.2), main = "b = 1, compensatory",  xlab = "Time step", ylab = "Population size [K]")
abline(h = 1, col = "black")
plot(gedaempft$sum..individuals..of.turtles[1:50] / 200, type = ltyg, ylim = c(0,2.2), main = "b = 2.5, damped oszillations",  xlab = "Time step", ylab = "")
abline(h = 1, col = "black")
plot(oszillierend$sum..individuals..of.turtles[1:50] / 200, type = ltyg, ylim = c(0,2.2), main = "b = 2.8, oscilating", xlab = "Time step", ylab = "Population size [K]")
abline(h = 1, col = "black")
plot(chaos$sum..individuals..of.turtles[1:50] / 200, type = ltyg, ylim = c(0,2.2), main = "b = 6, chaotic", xlab = "Time step", ylab = "")
abline(h = 1, col = "black")

#legend("bottomright", bg="white", inset=0.0, legend = c("species 1", "species2"), col = c("red", "blue"), merge = TRUE)

par <- oldpar

