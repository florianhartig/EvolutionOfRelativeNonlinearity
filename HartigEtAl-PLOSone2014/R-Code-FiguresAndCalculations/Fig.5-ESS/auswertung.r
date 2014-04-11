################################################################################
#
#
# http://florianhartig.wordpress.com/
################################################################################


rm(list=ls(all=TRUE)) 

require(fields)
require(graphics)
require(colorRamps)


wdir <- "/Users/Florian/Home/Projekte/Papers_published/14-HartigEtAl-PLOSone-SympatricEvolutionOfRNC/Code/Final/Figures/Fig.5-ESS"
setwd(wdir)


rescale_colors <- function(colorvector, breakpoint , pow1, pow2 ){
  
  len = length(colorvector)
  distr = breakpoint
  
  lower = round( breakpoint * len * (seq(0,1,length.out = (distr * 5000))^pow1)) 
  upper = round( breakpoint * len + (1-breakpoint) * len * (1- seq(1,0,length.out = ((1-distr) * 5000))^pow2 )) 

  out = colorvector[c(lower, upper)]
  return(out)
  
}


data <- read.csv("data/Coexistence ess-table.csv", skip = 6, )
data = data[order(data$X.run.number.),]
colnames(data)
critical <- read.csv("data/Coexistence ess-critical-table.csv", skip = 6, )
critical = critical[order(critical$X.run.number.),]
colnames(critical)

x1 <- unique(data$initial.rate.1) 
x2 <- unique(data$disturbance) 
c1 <- unique(critical$initial.rate.1) 
c2 <- unique(critical$disturbance) 

dx1 = x1[2] - x1[1]
dx2 = x2[2] - x2[1]

averageb <- (matrix(as.numeric(data$average.b), byrow=F, ncol=length(x1))) 
criticalb <- (matrix(as.numeric(critical$critical.b), byrow=F, ncol=length(c2)))

criticalb <- function(r,d){return(2-(2/(1+r*(d-1))))}
criticaltheoryvector <- criticalb(as.numeric(data$initial.rate.1), as.numeric(data$disturbance))
criticaltheoryvector[criticaltheoryvector < 0] = Inf
criticaltheoryvector[criticaltheoryvector > 25] = NaN
criticaltheorymatrix <- (matrix(criticaltheoryvector, byrow=F, ncol=length(x1))) 

#smoothened <- image.smooth(averageb, dx=dx1, dy = dx2, theta= .025 )


colorscheme = colorRampPalette(c("grey99", "yellow", "darkgreen", "black")) 
colorscheme2 = colorRampPalette(c("white","grey96", "black"))
colorscheme3 = colorRampPalette(c("lightgrey", "black"))
colorscheme4 = hcl(h = seq(210, 60, length = 4))
colorscheme5 = colorRampPalette(c("black", "black"))

colorscheme6 = colorRampPalette(c("darkred", "firebrick", "darkorange3", "orange", "gold", "yellow", "lemonchiffon", "lightyellow")) 

colorscheme7 = colorRampPalette(c("lemonchiffon", "yellow", " olivedrab1", "yellowgreen", "darkgreen", "black")) 
colorscheme8 = colorRampPalette(c("navyblue", "lightslateblue",  "grey96", "lightpink", "darkred"))


colorscheme9 = colorRampPalette(c("white", "lightpink", "palevioletred", "maroon4", "navyblue", "#000024")) 
#
#graphics.off()
#windows(width= 340, height = 310)
#oldpar <- par(mar = c(6,6,6,6), cex = 1.3, cex.lab = 1.4 ) # set graphic device
#

#graphics.off()
#windows(width= 600, height = 310)

#png("out.png", width= 800, height = 414)

pdf("fig5-raw.pdf", width = 20, height = 10)

par(mfrow= c(1,2), cex.lab = 1.2)



#globalmin = min(c(min(averageb),min(criticalb))) 
#globalmax = max(c(max(averageb))) 

col= colorscheme(100)
basecol = matlab.like(100000)
basecol = rescale_colors(basecol, 0.41,3,3)

col = c(colorscheme5(100), colorscheme9(3000)[3000:1])
#col = rescale_colors(basecol, 0.75, 3, 3)
logfactor = log(10)
par(mar = c(6,6,6,6), cex = 1.4)
image.plot(x1, x2, averageb,  col = col,  main = "(a) ESS of individual-based model", font.main=1, xlab = "Intrinsic growth rate, r", ylab = "Disturbance mortality, d")      #
#contour(x1, x2, smoothened$z, levels = c(3,5,10), add = T)  
par(mar = c(6,6,6,6), cex = 1.4)
#image.plot(c1, c2, criticalb, col = col, zlim = zlimits,  main = "Critical b-strategy\n(Classical MSS model)", xlab = "Growth rate r", ylab = "Disturbance d")      #
image.plot(x1, x2, averageb-criticaltheorymatrix, col = colorscheme8(100)[100:1], zlim = c(-10, 10), main = "(b) Deviation of ESS from critical b-value", font.main=1, xlab = "Intrinsic growth rate, r", ylab = "Disturbance mortality, d")      #

                     
#dev.copy(png,'out.png')
dev.off()


