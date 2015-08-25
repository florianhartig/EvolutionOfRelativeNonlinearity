#Population function---------

msslvsepopfun <- function(parmat, start, time ){
  
  #Setting required Vectors and Parameters
  N1<- numeric(time)
  N2 <- numeric(time)
  N1[1] <- start[1]
  N2[1] <- start[2]
  
  fn1 <- numeric(time)
  fn2 <- numeric(time)
  
  dn1 <- numeric(time)
  dn2 <- numeric(time)
  
  r.capita1 <- numeric(time)
  r.capita2 <- numeric(time)
  
  buffer1 <- numeric(time)
  buffer2 <-numeric(time)
  buffer1[1] <- 250000
  buffer2[1] <- 250000
  
  
  
  #Getting Parameter Values from the Parametermatrix   
  r1=parmat[1,1]
  d1=parmat[2,1]
  alpha12 = parmat[3,1]
  k1=parmat[4,1]
  b1=parmat[5,1]
  r2=parmat[1,2]
  d2=parmat[2,2]
  alpha21= parmat[3,2]
  k2=parmat[4,2]
  b2=parmat[5,2]
  
  #Drawing environmental stochasticity
  env <- rlnorm(time,meanlog =0, sdlog =parmat[6,1]) 
  
  #Running the for loop for discrete population growth
  
  for (i in c(1:(time-1))){
    
    #If the storage Effect is turned on the environmental responses of the competing #species are opposing 
    #("different environmental response")
    if(parmat[7,1]==1){env1 <- env[i]
                       env2 <- 1/env[i]}
    else{env1 <- env[i]
         env2 <- env[i]}
    
    a1 <- (r1)  / (1 + (r1 -1) * ((N1[i] +  alpha12*N2[i])/k1*env1)^b1) 
    a2 <- (r2)  / (1 + (r2 -1) * ((N2[i]  + alpha21*N1[i])/k2*env2)^b2) 
    
    #Just in case the reproduction ratio gets negative it is set to zero here. This makes sense because later it is multiplied with N(t) and negative populations do not make sense. In this case the population is extinct. Also it would create NA's when drawing from the poisson distribution.
    a1[a1<0]<-0
    a2[a2<0]<-0
    #If the storage Effect is turned on, buffered population growth is modeled here. If #not a1 and a2 are
    #the reproduction ratios and can just be multiplied with N(t)
    if(parmat[7,1]==1)
    {
      
      
      
      
      #Filling the buffer (for example producing seeds for a seed bank)
      buffer1[i] <- buffer1[i] + a1*N1[i]*1000
      buffer2[i] <- buffer2[i] + a2*N2[i]*1000
      
      #The rate how many seeds start to germinate and build up the new generation is #environment dependent
      N1[i+1] <-  0.00005*env1* buffer1[i]
      N2[i+1] <-  0.00005*env2* buffer2[i]
      
      #The germinated seeds as well as the destroyed seeds (factor 0.2) are subtracted #from the current buffer
      buffer1[i+1] <- buffer1[i] - N1[i+1]
      buffer1[i+1] <- buffer1[i] - 0.2*buffer1[i]
      
      buffer2[i+1] <- buffer2[i] - N2[i+1]
      buffer2[i+1] <- buffer2[i] - 0.2*buffer2[i]
    }
    else{
      #When storage effect is not ivolved the adults who died are subtracted here (in the #storage effect this
      # has already been taken into account with the factor for destroyed seeds)
      N1[i+1] <- N1[i]*(1-d1) *a1
      N2[i+1] <- N2[i]*(1-d2) *a2
    }
    #For demographic stochasticity the resulting population density for the next time step is drawn from a poisson       
    #distribution
    
    
    N1[i+1] <- rpois(1,N1[i+1])
    N2[i+1] <- rpois(1,N2[i+1])
    
    
    
    #storing the reproduction ratio
    fn1[i] <- N1[i+1]/N1[i]
    fn2[i] <- N2[i+1]/N2[i]
    
  }
  #Calculating delta N
  dn1 <- N1-c(0,N1[-time]) 
  dn2 <- N2-c(0,N2[-time])
  dn1[1] <- NA
  dn2[1] <- NA
  
  #Calculating the per capita growth rate
  for(k in c(1:time-1)){
    
    r.capita1[k] <- dn1[k+1]/N1[k]
    r.capita2[k] <- dn2[k+1]/N2[k]
    
  }
  popmat <- cbind(N1,N2,fn1,fn2,dn1,dn2,r.capita1,r.capita2,buffer1,buffer2,env)
  
  return(popmat)
}

#Parmats-------
parmatneutral <- matrix(ncol = 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatneutral[,1] <- c(5,0.05,1,500,1,0.15,0)
parmatneutral[,2] <- c(5,0.05,1,500,1,0.15,0)

parmatrnc <- matrix(ncol = 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatrnc[,1] <- c(5,0.05,1,500,5,0.15,0)
parmatrnc[,2] <- c(5,0.05,1,500,1,0.15,0)

parmatlv <- matrix(ncol = 2, nrow = 7)
#F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
#Parametersequence: r,d,alpha,k,b,env,se
parmatlv[,1] <- c(5,0.05,0.6,500,1,0.15,0)
parmatlv[,2] <- c(5,0.05,0.9,500,1,0.15,0)

parmatse <- matrix(ncol = 2, nrow = 7)
#F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
#Parametersequence: r,d,alpha,k,b,env,se
parmatse[,1] <- c(5,0.05,1,500,1,1.2,1)
parmatse[,2] <- c(5,0.05,1,500,1,1.2,1)

parmatcomb1 <- matrix(ncol= 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatcomb1[,1] <- c(5,0.05,1,500,5,1.2,1)
parmatcomb1[,2] <- c(5,0.05,1,500,1,1.2,1)

parmatcomb2 <- matrix(ncol= 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatcomb2[,1] <- c(5,0.05,0.6,500,5,0.15,0)
parmatcomb2[,2] <- c(5,0.05,0.9,500,1,0.15,0)

parmatcomb3 <- matrix(ncol= 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatcomb3[,1] <- c(5,0.05,0.9,500,1,1.2,1)
parmatcomb3[,2] <- c(5,0.05,0.6,500,1,1.2,1)

parmatall <- matrix(ncol= 2, nrow = 7)
#Parametersequence: r,d,alpha,k,b,env,se
parmatall[,1] <- c(5,0.05,0.6,500,5,1.2,1)
parmatall[,2] <- c(5,0.05,0.9,500,1,1.2,1)

#Fluctuation imbalance-----
In <- numeric(1000) 
Ir <- numeric(1000)
Il <- numeric(1000)
Is <- numeric(1000)
Ic1 <-numeric(1000)
Ic2 <-numeric(1000)


for(i in 1:1000){
  #Neutral
  pop <- msslvsepopfun(parmat=parmatneutral,start=c(200,200), time = 500)
  In[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))  
  #RNC
  pop <- msslvsepopfun(parmat=parmatrnc,start=c(200,200), time = 500)
  Ir[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))  
  #Niche
  pop <- msslvsepopfun(parmat=parmatlv,start=c(200,200), time = 500)
  Il[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))
  #Storage Effect
  pop <- msslvsepopfun(parmat=parmatse,start=c(200,200), time = 500)
  Is[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))  
  #Combination 1
  pop <- msslvsepopfun(parmat=parmatcomb1,start=c(200,200), time = 500)
  Ic1[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))  
  #Combination 2
  pop <- msslvsepopfun(parmat=parmatcomb2,start=c(200,200), time = 500)
  Ic2[i] <-  length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))*length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))  
}

#Fluctuation difference--------
Dn <- numeric(1000)
Dr<- numeric(1000)
Dl<- numeric(1000)
Ds<- numeric(1000)
Dc1<- numeric(1000)
Dc2<- numeric(1000)

for(i in c(1:1000)){
  #Neutral
  pop <- msslvsepopfun(parmat=parmatneutral,start=c(200,200), time = 500)
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Dn[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  #RNC
  pop <- msslvsepopfun(parmat=parmatrnc,start=c(200,200), time = 500)
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Dr[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  #Niche Partitioning
  pop <- msslvsepopfun(parmat=parmatlv,start=c(200,200), time = 500)
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Dl[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  #Storage Effect
  pop <- msslvsepopfun(parmat=parmatse,start=c(200,200), time = 500 )
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Ds[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  #Combination 1
  pop <- msslvsepopfun(parmat=parmatcomb1, start=c(200,200), time = 500 )
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Dc1[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  #Combination 2
  pop  <- msslvsepopfun(parmat=parmatcomb2, start=c(200,200), time = 500)
  pop[pop[,1]==0]<-NA
  pop[pop[,2]==0]<-NA
  Dc2[i] <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  
}

#Correlation-----
rn <- numeric(1000)
rl <- numeric(1000)


for(i in c(1:1000)){
  
  
  neutral <- msslvsepopfun(parmat=parmatneutral,start=c(200,200), time = 500)
  lv <- msslvsepopfun(parmat=parmatlv,start=c(200,200), time = 500)
  
  rn[i] <- abs(cor(neutral[,1],neutral[,2]))
  rl[i] <- abs(cor(lv[,1],lv[,2]))
  
}

#Mean analysis------
Ms <- numeric(1000)
Mc1 <- numeric(1000) 

for(i in c(1:1000)){
  
  
  se <- msslvsepopfun(parmat=parmatse,start=c(200,200), time = 500)
  comb1 <- msslvsepopfun(parmat=parmatcomb1, start=c(200,200), time = 500)
  
  Ms[i] <- abs(1- mean(se[,1])/mean(se[,2]))
  Mc1[i] <- abs(1- mean(comb1[,1])/mean(comb1[,2]))
  
}

#Fluctuation frequency---------
Fr <- numeric(1000)
Fc2 <- numeric(1000)

for(i in 1:1000){
  
  rnc <- msslvsepopfun(parmat=parmatrnc,start=c(200,200),time = 500)
  comb2 <- msslvsepopfun(parmat=parmatcomb2, start=c(200,200), time = 500)
  
  
  Fr[i] <-length(which(abs(rnc[,5])<mean(rnc[,1])/2))/length(rnc)
  
  
  Fc2[i] <-length(which(abs(comb2[,5])<mean(comb2[,1])/2))/length(rnc)
  
  
}

#Validation of threshold values-----
binom.test(length(which(Is>4.5))+ length(which(Ic1 > 4.5)),2000)
binom.test(length(which(In<4.5))+ length(which(Ir<4.5))+length(which(Il<4.5))+length(which(Ic2<4.5)),4000)

binom.test(length(which(Dr>1))+length(which(Dc2>1)), 2000)
binom.test(length(which(Dn<1))+length(which(Dl<1))+length(which(Ds<1))+length(which(Dc1<1)), 4000)

binom.test(length(which(rn>0.5)), 1000)
binom.test(length(which(rl<0.5)), 1000)

binom.test(length(which(Mc1>0.9)),1000)
binom.test(length(which(Ms <0.9)),1000)

binom.test(length(which(Fc2<0.03)),1000)
binom.test(length(which(Fr >0.03)),1000)

#Exclusion Approach-----
exclusion <- function(TIMES,K=F,B=F,Alpha=F){
  
  #For later comparison the information which mechanism is drawn gets stored in "Truth". In the vector Guess the result of the analysation is stored.
  Truth <- character(TIMES)
  Guess <- character(TIMES)
  mechanisms <- c("neutral", "rnc", "lv", "se", "comb1" ,"comb2")
  storage= character(1)
  
  for(i in c(1:TIMES)){
    
    #First the parameter values are drawn from their specific interval if they are to be alternated.
    if(K==T){capacity <-  sample(seq(from=500, to = 10000, by = 100),size = 1)}
    else{capacity <- 500}
    
    if(B==T){
      b1 <- sample(c(seq(from=0.4, to = 1, by = 0.1),seq(from=3, to = 7, by =0.2)),size = 1)
      if(b1>2.9){ b2 <- sample(seq(from=0.4, to = 1, by =0.1), size =1)}
      else{b2 <- sample(seq(from = 3, to =7, by = 0.2),size=1)}
    }    
    else{b1 <- 5
         b2 <- 1}
    
    if(Alpha==T){alpha1 <- sample(seq(from=0.1, to = 0.9, by =0.02),size=1) 
                 alpha2 <- sample(seq(from=0.1, to = 0.9, by =0.02),size=1)}
    else{alpha1 <- 0.6
         alpha2 <- 0.9}
    
    #Then the values are stored in the parmats
    
    parmatneutral <- matrix(ncol = 2, nrow = 7)
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatneutral[,1] <- c(5,0.05,1,capacity,1,0.15,0)
    parmatneutral[,2] <- c(5,0.05,1,capacity,1,0.15,0)
    
    parmatrnc <- matrix(ncol = 2, nrow = 7)
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatrnc[,1] <- c(5,0.05,1,capacity,b1,0.15,0)
    parmatrnc[,2] <- c(5,0.05,1,capacity,b2,0.15,0)
    
    parmatlv <- matrix(ncol = 2, nrow = 7)
    #F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatlv[,1] <- c(5,0.05,alpha1,capacity,1,0.15,0)
    parmatlv[,2] <- c(5,0.05,alpha2,capacity,1,0.15,0)
    
    parmatse <- matrix(ncol = 2, nrow = 7)
    #F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatse[,1] <- c(5,0.05,1,capacity,1,1.2,1)
    parmatse[,2] <- c(5,0.05,1,capacity,1,1.2,1)
    
    parmatcomb1 <- matrix(ncol= 2, nrow = 7)
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatcomb1[,1] <- c(5,0.05,1,capacity,b1,1.2,1)
    parmatcomb1[,2] <- c(5,0.05,1,capacity,b2,1.2,1)
    
    parmatcomb2 <- matrix(ncol= 2, nrow = 7)
    #Parameterreihenfolge: r,d,alpha,k,b,env,se
    parmatcomb2[,1] <- c(5,0.05,alpha1,capacity,b1,0.15,0)
    parmatcomb2[,2] <- c(5,0.05,alpha2,capacity,b2,0.15,0)
    
    #Here the mechanism for the current simulation is drawn and simulated
    Truth[i] <- sample(mechanisms,size=1)
    if(Truth[i] =="neutral"){parmat = parmatneutral}
    if(Truth[i] == "rnc")   {parmat = parmatrnc}
    if(Truth[i] == "lv")    {parmat = parmatlv}
    if(Truth[i] == "se")    {parmat = parmatse}
    if(Truth[i] == "comb1") {parmat = parmatcomb1}
    if(Truth[i] == "comb2") {parmat = parmatcomb2}
    
    pop <-msslvsepopfun(parmat=parmat,start=c(200,200), time = 500)
    
    #The 5 assessment tools are calculated    
    i1 <- length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))
    i2 <- length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))
    I  <- i1*i2  
    
    r <- cor(pop[,1],pop[,2])
    k <- mean(pop[,1])+ mean(pop[,2])
    rcrit <- k^(1/5.2877125)/3.0500907 - 0.9590107
    
    f1 <- length(which(abs(pop[,5])<mean(pop[,1])/2))/length(pop)    
    f2 <- length(which(abs(pop[,6])<mean(pop[,2])/2))/length(pop)    
    
    M <- abs(1- mean(pop[,1])/mean(pop[,2]))
    
    #Because in the calculation of the fluctuation difference is a division with the population density involved, we need to make sure they are not zero (e.g. if one species goes extinct). Otherwise the result would be an infinite value.  
    pop[,1][pop[,1] == 0] <- NA
    pop[,2][pop[,2] == 0] <- NA
    D <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
    
    #Now the exclusion of the mechanisms begins(Following the schema). 
    if(I>4.5){Guess[i] <- "se"}
    else{    
      if(D > 1){ Guess[i] <- "rnc"}
      else{
        if(K!=T){if(abs(r) > 0.5){Guess[i] <- "neutral"}else{Guess[i] <- "lv"}}
        else{if(r < rcrit){Guess[i] <- "neutral"}else{Guess[i] <- "lv"}}
      }
    }
    
    
    #After the involving mechanisms are identified it is possible to determine Combinations.
    if(Guess[i] == "se"){if(M > 0.9){Guess[i] <- "comb1"}}
    if(Guess[i] == "rnc"){if(f1 < 0.03){Guess[i] <- "comb2"}
                          if(f2 < 0.03){Guess[i] <- "comb2"}}
    
  }
  result<- cbind(Truth,Guess)
  return(result)
}

#Validating exclusion approach------
x <- exclusion(1000)
binom.test(length(which(x[,1]==x[,2])),1000)

x <- exclusion(1000,K=T)
binom.test(length(which(x[,1]==x[,2])),1000)

x <-exclusion(1000,K=T,Alpha=T)
binom.test(length(which(x[,1]==x[,2])),1000)

x<-exclusion(1000,K=T,B=T)
binom.test(length(which(x[,1]==x[,2])),1000)

x <- exclusion(1000,B=T,Alpha=T)
binom.test(length(which(x[,1]==x[,2])),1000)

x<-exclusion(1000,K=T,B=T,Alpha=T)
binom.test(length(which(x[,1]==x[,2])),1000)

#Classifier approach-------

classifierfun <- function(mechanism, K=F,Alpha=F,B=F){
  
  #First the parameter values are drawn from their specific interval if they are to be alternated.
  if(K==T){capacity <-  sample(seq(from=500, to = 10000, by = 100),size = 1)}
  else{capacity <- 500}
  
  if(B==T){
    b1 <- sample(c(seq(from=0.4, to = 1, by = 0.1),seq(from=3, to = 7, by =0.2)),size = 1)
    if(b1>2.9){ b2 <- sample(seq(from=0.4, to = 1, by =0.1), size =1)}
    else{b2 <- sample(seq(from = 3, to =7, by = 0.2),size=1)}
  }    
  else{b1 <- 5
       b2 <- 1}
  
  if(Alpha==T){alpha1 <- sample(seq(from=0.1, to = 0.9, by =0.02),size=1) 
               alpha2 <- sample(seq(from=0.1, to = 0.9, by =0.02),size=1)}
  else{alpha1 <- 0.6
       alpha2 <- 0.9}
  
  #Then the values are stored in the parmats
  
  parmatneutral <- matrix(ncol = 2, nrow = 7)
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatneutral[,1] <- c(5,0.05,1,capacity,1,0.15,0)
  parmatneutral[,2] <- c(5,0.05,1,capacity,1,0.15,0)
  
  parmatrnc <- matrix(ncol = 2, nrow = 7)
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatrnc[,1] <- c(5,0.05,1,capacity,b1,0.15,0)
  parmatrnc[,2] <- c(5,0.05,1,capacity,b2,0.15,0)
  
  parmatlv <- matrix(ncol = 2, nrow = 7)
  #F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatlv[,1] <- c(5,0.05,alpha1,capacity,1,0.15,0)
  parmatlv[,2] <- c(5,0.05,alpha2,capacity,1,0.15,0)
  
  parmatse <- matrix(ncol = 2, nrow = 7)
  #F¸r Stable : K1 > K2 * alpha12 und K2 > K1*alpha21
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatse[,1] <- c(5,0.05,1,capacity,1,1.2,1)
  parmatse[,2] <- c(5,0.05,1,capacity,1,1.2,1)
  
  parmatcomb1 <- matrix(ncol= 2, nrow = 7)
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatcomb1[,1] <- c(5,0.05,1,capacity,b1,1.2,1)
  parmatcomb1[,2] <- c(5,0.05,1,capacity,b2,1.2,1)
  
  parmatcomb2 <- matrix(ncol= 2, nrow = 7)
  #Parameterreihenfolge: r,d,alpha,k,b,env,se
  parmatcomb2[,1] <- c(5,0.05,alpha1,capacity,b1,0.15,0)
  parmatcomb2[,2] <- c(5,0.05,alpha2,capacity,b2,0.15,0)
  
  
  if(mechanism =="neutral") {parmat = parmatneutral}
  if(mechanism == "rnc")    {parmat = parmatrnc}
  if(mechanism == "lv")     {parmat = parmatlv}
  if(mechanism == "se")     {parmat = parmatse}
  if(mechanism == "comb1")  {parmat = parmatcomb1}
  if(mechanism == "comb2")  {parmat = parmatcomb2}
  
  #Every mechanism gets its own class  
  if(mechanism =="neutral") {class <- 1}
  if(mechanism == "rnc")    {class <- 2}
  if(mechanism == "lv")     {class <- 3}
  if(mechanism == "se")     {class <- 4}
  if(mechanism == "comb1")  {class <- 5}
  if(mechanism == "comb2")  {class <- 6}
  
  a <- numeric(7)
  pop <-msslvsepopfun(parmat=parmat,start=c(200,200), time = 500)
  
  #The 5 assessment tools are calculated    
  i1 <- length(which(pop[,1]<mean(pop[,1])))/length(which(pop[,1]>mean(pop[,1])))
  i2 <- length(which(pop[,2]<mean(pop[,2])))/length(which(pop[,2]>mean(pop[,2])))
  I  <- i1*i2  
  
  r <- cor(pop[,1],pop[,2])
  
  f1 <- length(which(abs(pop[,5])<mean(pop[,1])/2))/length(pop)    
  f2 <- length(which(abs(pop[,6])<mean(pop[,2])/2))/length(pop)    
  
  M <- abs(1- mean(pop[,1])/mean(pop[,2]))
  
  #Because in the calculation of the fluctuation difference is a division with the population density involved, we need to make sure they are not zero (e.g. if one species goes extinct). Otherwise the result would be an infinite value.  
  pop[,1][pop[,1] == 0] <- NA
  pop[,2][pop[,2] == 0] <- NA
  D <- abs(mean(abs(pop[,5])/(pop[,1]),na.rm=T) - mean(abs(pop[,6])/pop[,2],na.rm=T))
  
  
  
  a[1] <- r
  a[2] <- D
  a[3] <- I
  a[4] <- f1
  a[5] <- f2
  a[6] <- M
  a[7] <- class
  
  return(a)
}


classmatrix <- matrix(ncol=7,nrow=600)

for(i in 1:100){
  classmatrix[i,] <- classifierfun(mechanism="neutral")
}
for(i in 1:100){
  classmatrix[100+i,] <- classifierfun(mechanism="rnc")
}
for(i in 1:100){
  classmatrix[200+i,] <- classifierfun(mechanism="lv")
}
for(i in 1:100){
  classmatrix[300+i,] <- classifierfun(mechanism="se")
}
for(i in 1:100){
  classmatrix[400+i,] <- classifierfun(mechanism="comb1")
}
for(i in 1:100){
  classmatrix[500+i,] <- classifierfun(mechanism="comb2")
}


install.packages("caret")
require("caret")

TrainData <- classmatrix[,1:6]
TrainClasses <- classmatrix[,7]

cforestfit <- train(TrainData, TrainClasses,
                      method = "cforest")

#Validating Classifier approach---------

mat <- matrix(ncol = 7, nrow=1000)
for(i in 1:1000){
  mat[i,]<- classifierfun(mechanism=sample(c("neutral", "rnc", "lv", "se", "comb1", "comb2"),size=1))
  }

result <- round(predict(cforestfit, newdata = mat[,1:6]))
x <- cbind(result,mat[,7])

binom.test(length(which(x[,1]==x[,2])),1000)

