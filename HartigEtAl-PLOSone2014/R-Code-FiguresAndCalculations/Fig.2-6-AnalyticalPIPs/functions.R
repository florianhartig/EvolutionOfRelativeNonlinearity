################################################################################
#
# Florian Hartig                             http://florianhartig.wordpress.com/
################################################################################

library(colorRamps)
require(compiler)


mss <- function(N, r=5, b=2, K=1, d=0.05, procerror = F, max = 0.05){
  if (procerror==T){
    return(r * (1-d) / (1 + (r-1) * (N/K)^b) * (1+runif(1,max=max)-max/2)  )
  }
  else{
    return(r * (1-d) / (1 + (r-1) * (N/K)^b) )
  }
}

mss <- cmpfun(mss)




mssmod <- function(N, r=5, b=2, K=1, d=0.05, scaling = 4, lower=0.8, upper=5){
  
  middle = upper # middle is not used any more, so always upper = middle
  
  if (b<=lower || b>=upper){
    return(r * (1-d) / (1 + (r-1) * (N/K)^b) )
  }
  
  else if (b<middle){
    b=ifelse(N < 1,
             lower + (b-lower)^scaling / (middle-lower)^(scaling-1),
             lower + (b-lower)^(1/scaling) / (middle-lower)^(1/scaling-1)
    )
    return(r * (1-d) / (1 + (r-1) * (N/K)^b) )
  }
  #   else if (b<upper){
  #     b=ifelse(N > 1,
  #              middle + (b-middle)^scaling / (upper-middle)^(scaling-1),
  #              middle + (b-middle)^(1/scaling) / (upper-middle)^(1/scaling-1)
  #     )
  #     return(r * (1-d) / (1 + (r-1) * (N/K)^b) )
  #   }
}

smod <- cmpfun(mssmod)




ricker <- function(N, r=0.5, b=2, K=1, d=0.05){
  
  return((1-d) * exp( r * (1-(N/K)^b)))

}

ricker <- cmpfun(ricker)


hassell <- function(N, r=40, b=2, K=1, d=0.05){

  return((1-d) * r  / (1 + (r^( 1/b) - 1) * (N / K))^b)
  
}

hassell <- cmpfun(hassell)



create_resident_dynamics <- function(popfunction, startvalue=0.9, n=500){
  out = rep(NA,(n+101))
  out[1]=startvalue
  for (i in 1:(n+100)){
    out[i+1] = out[i] * popfunction(out[i])
  }
  return(out[102:(101+n)])
}




plot_popdyn <- function(testpop){
  par(mfrow = c(2,1))
  plot(testpop, type="b")
  hist(testpop, breaks = 50)
  abline(v=mean(testpop), col = "red")
}



basecol = matlab.like(1000000)


rescale_colors <- function(colorvector, breakpoint , pow1, pow2 ){
  
  #breakpoint = ( abs(minmax[1]) - abs(minmax[2] - minmax[1]) ) / abs(minmax[2] - minmax[1])
  
  len = length(colorvector)
  distr = breakpoint
  
  lower = round( breakpoint * len * (seq(0,1,length.out = (distr * 5000))^pow1)) 
  upper = round( breakpoint * len + (1-breakpoint) * len * (1- seq(1,0,length.out = ((1-distr) * 5000))^pow2 )) 
  
  out = colorvector[c(lower, upper)]
  return(out)
  
}




create_PIP <- function(popfun, residentcycles = 500, spacing = 200){
  #bval = exp(seq(log(0.1),log(20), length.out=spacing))
  bval = exp(seq(log(0.1),log(20), length.out=spacing))
  invasionfitness = matrix(NA, nrow = spacing, ncol = spacing)
  for (i in 1:length(bval)){
    newfunc <- function(n)popfun(n,b=bval[i])
    
    residents <- create_resident_dynamics(newfunc,n=residentcycles)
    #debuggingn
    #plot_popdyn(residents)
    for (j in 1:length(bval)){
      newfunc <- function(n)popfun(n,b=bval[j])
      invasionfitness[i,j] <- sum(log(sapply(residents, newfunc))) / residentcycles
    }
  }
  return(list(bval, invasionfitness))
}

