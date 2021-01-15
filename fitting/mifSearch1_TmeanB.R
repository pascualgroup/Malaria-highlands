rm(list=ls())

require(pomp)
require(plyr)

dat <- read.csv("DataCovar_1975-2000.csv",header=T)	# read in data file with the malaria data and covariates

x <- as.numeric(Sys.getenv("PBS_ARRAYID"))
seqP <- seq(from = 1, to = 50000, by = 100)
now.num <- seqP[x]


y <- read.csv("parameters_TmeanB.csv")
param <- as.numeric(y[now.num,])
param.names <- colnames(y)
names(param) <- param.names

source("poObject_TmeanB.R")



for(cont in now.num:(now.num+99)){
  for(i in 1:3){

	seed <- ceiling(runif(1,min=1,max=2^30))
	set.seed(seed)
	param <- as.numeric(y[cont,])
        names(param) <- param.names
        cat(cont, i, "\n")


tryCatch(mif2(
      po,
      Nmif = 50,              #change from 10 to 50
      start = param,
      Np = 1000, #change from 1000 to 15000
      cooling.type="hyperbolic",
      cooling.fraction.50 = 0.5,
      rw.sd = rw.sd(muEI = 0.03, muIQ = 0.03, muIS = 0.03, muQS = 0.03,  sigPRO = 0.03, sigOBS = 0.03,
                    tau = 0.03, rho = 0.03, betaOUT = 0.03, bT4 = 0.03, bT6 = 0.03, b1 = 0.03, b2 = 0.03, 
                    b3 = 0.03, b4 = 0.03, b5 = 0.03, b6 = 0.03, q0 = 0.03,
                    S_0 = 0.03, E_0 = 0.03, I_0 = 0.03, Q_0 = 0.03 , K_0 = 0.03, F_0 = 0.03),
      transform=TRUE),  error = function(e) e) -> mifout

    if(length(coef(mifout)) > 0){
      loglik.mif <- replicate(n=5,logLik(pfilter(po,
                    params=coef(mifout),Np=1000,max.fail=500)))
      bl <- logmeanexp(loglik.mif,se=TRUE)
      loglik.mif.est <- bl[1]
      loglik.mif.se <- bl[2]
      cat(cont,loglik.mif.est,"\n")

 ### SAVE OUTPUT
      if(is.finite(loglik.mif.est)){
        par.out <- coef(mifout)
        names(par.out) <- param.names
        if(file.exists("mifOutput_TmeanB.csv")){
          write.table(t(as.matrix(c(cont,seed,par.out,loglik.mif.est,loglik.mif.se))), 
                      "mifOutput_TmeanB.csv", append = T, 
                      col.names=F, row.names=F, sep = ",")
        } else{ 
          write.table(t(as.matrix(c(cont,seed,par.out,loglik.mif.est,loglik.mif.se))), 
                      "mifOutput_TmeanB.csv", append = T, col.names=c("run","seed", 
                                                                    param.names, "loglik", "loglik.se"), row.names=F, sep = ",") }
      }
    }
  }
}


 
