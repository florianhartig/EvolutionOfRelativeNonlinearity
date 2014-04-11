################################################################################
# Script for creating bifurkation plot in fig.1
#
# http://florianhartig.wordpress.com/
################################################################################

rm(list=ls(all=TRUE))

#install.packages("tgp")
require(fields)

      
wdir <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.1-ResponseAndBifurcation/"
setwd(wdir)           
                    
population <- read.csv(paste(wdir,"simulationOutput/Coexistence complexity-table.csv", sep=""), skip = 6, )

selection = population[order(population$X.run.number.),] 
selection = selection[selection$ticks > 1000,]
bvalues = unique(selection$initial.strategy.1)[100:200]


plot.new()

plot(0,0,  xaxt="n", xlim = c(log10(min(bvalues)), log10(max(bvalues))), ylim = c(0,2.5), 
		col = "white", xlab = "b-strategy", ylab = "Population size [K]")
axis(side = 1, at = log10(c(seq(from = 1, to = 7.38, by = 0.2), seq(from = 2, to = 9, by = 1), seq(from = 10 ,to = 50, by = 10) )), labels = F)
axis(side = 1, at = log10(c(1, 2, 4, 7)), labels = c("1", "2", "4", "7"), lwd.ticks = 2, tcl = - 0.5 )



chooser = T
bifvalue = 2.53
for (bval in bvalues){
  values = unique(round(selection[selection$initial.strategy.1 == bval, ]$sum..individuals..of.species, digits =3)) /200
  points(rep(log10(bval), length(values)), values, cex = 0.5)
}

abline(v= log10(bifvalue))




