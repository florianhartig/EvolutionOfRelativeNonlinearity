################################################################################
#
#
# http://florianhartig.wordpress.com/
################################################################################


library(fields)

################################################################################
# read in all data 

path <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.4-EvolutionCollapse"

setwd(path)

timesteps = 200
totaltime = 150000
lines = totaltime / timesteps + 1

bval = seq(0,7, length = 140)
timeval = seq(0,150000, by=200)  / 1000

bdouble = read.table("data/branching-double.dat", header = T, na.strings = "0")
blow = read.table("data/branching-low.dat", header = T, na.strings = "0")
bup = read.table("data/branching-up.dat", header = T, na.strings = "0")

mbdouble = as.matrix(bdouble[1:lines,])
mbup = as.matrix(bup[1:lines,])
mblow = as.matrix(blow[1:lines,]    )

################################################################################

pdf("fig4-raw.pdf", width = 20, height = 10)


par(mfrow = c(1,2))
par(cex.lab = 1.2)



pal = colorRampPalette(c( "orange", "yellow", "white"), bias = 1)
image(timeval, bval, mbdouble , col = pal(100), main = "Coexisting runs", font.main=1, xlab = "10^3 generations", ylab = "Density compensation strategy b")
abline(v = 30)

par(cex.lab = 1.2)
pal = colorRampPalette(c("orange", "yellow", "white"), bias = 1)
image(timeval, bval, mblow, col = pal(100), add = F, main = "Independent runs", xlab = "10^3 generations")

par(cex.lab = 1.2)
pal = colorRampPalette(c("darkolivegreen1", "yellow", "white"), bias = 1)
image(timeval, bval, mbup, col = pal(100), add = T)
abline(v = 30)


dev.off()

