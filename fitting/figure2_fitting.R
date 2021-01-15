################################
require(plyr)
require(pomp)
require(ggplot2)
require(reshape2)

dat <- read.csv("DataCovar_1975-2000.csv",header=T)
Tmean <- read.csv("mifOutputProf_rho.csv")
Tmean <- subset(Tmean, !is.na(loglik))
Tmean <- subset(Tmean, !is.na(loglik.se))
Tmean <- arrange(Tmean, -loglik)
Tmean <- Tmean[,2:27]

now.num <- 1
param <- as.numeric(Tmean[now.num,])
param.names <- colnames(Tmean)
names(param) <- param.names

source("poObject_TmeanB.R")
datL <- read.csv("DataCovar_1975-2007.csv",header=T)	# read in data file with the malaria data and covariates
datO <- subset(dat, select=c("time", "Pf"))

simT <- NULL
for(j in 1:1000){
  sim.test <- simulate(po,params=param,nsim=1)
  simT <- cbind(simT, data.frame(sim.test)[,'Pf'])
}
sim.test <- simulate(po,params=param,nsim=1)
sim.test <- data.frame(sim.test)

dat$mean <- apply(simT, 1, function(x){mean(x, na.rm = T)})
dat$median <- apply(simT, 1, function(x){median(x, na.rm = T)})
dat$qtLow <- apply(simT, 1, function(x){quantile(x,probs = 0.1, na.rm = T) })
dat$qtHig <- apply(simT, 1, function(x){quantile(x,probs = 0.9, na.rm = T) })

pfout <- pfilter(po,params=param,Np=1000,save.states=T)

S0  <- pfout@saved.states[[length(pfout@saved.states)]]['S',]
E0  <- pfout@saved.states[[length(pfout@saved.states)]]['E',]
I0  <- pfout@saved.states[[length(pfout@saved.states)]]['I',]
Q0  <- pfout@saved.states[[length(pfout@saved.states)]]['Q',]
K0  <- pfout@saved.states[[length(pfout@saved.states)]]['K',]
F0  <- pfout@saved.states[[length(pfout@saved.states)]]['F',]

totpop  <- S0[1] + E0[1] + I0[1] + Q0[1] 
S0 <- S0/totpop
E0 <- E0/totpop
I0 <- I0/totpop
Q0 <- Q0/totpop

datfut <- subset(datL, select = time)
datfut <- subset(datfut, time > 2000)
datfut$Pf <- NA

po2 <- pomp(
  data=datfut,     
  times='time',
  t0=2000,
  covar=subset(datL,select=c('time','tcovn2','Tm1t','Tm2t','pop','dpopdt',sprintf('season%d',1:6))),
  rprocess=euler.sim(Csnippet(simul),delta.t=2/365),
  rmeasure = Csnippet(rmeas),
  dmeasure = Csnippet(dmeas),
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  initializer=initlz,
  tcovar = 'time',
  zeronames = c('cases','W','err'),
  statenames = c('cases','err','S','E','I','Q','K','F','W'),
  paramnames = names(param)
)

simT <- NULL
for(i in 1:1000){
  param[c('S_0','E_0','I_0','Q_0','K_0','F_0')] <- c(S0[i],E0[i],I0[i],Q0[i],K0[i],F0[i])
  sim.test <- simulate(po2,params=param,nsim=1)
  simT <- cbind(simT, data.frame(sim.test)[,'Pf'])
}
sim.test <- simulate(po2,params=param,nsim=1)
sim.test <- data.frame(sim.test)

datfut$mean <- apply(simT, 1, function(x){mean(x, na.rm = T)})
datfut$median <- apply(simT, 1, function(x){median(x, na.rm = T)})
datfut$qtLow <- apply(simT, 1, function(x){quantile(x,probs = 0.1, na.rm = T) })
datfut$qtHig <- apply(simT, 1, function(x){quantile(x,probs = 0.9, na.rm = T) })

dat <- subset(dat, select = colnames(datfut))
dat <- rbind(dat, datfut[1,])

########## Plot ########## 

p1 <- ggplot() + theme_classic() +  xlab("Year") + ylab("Reported cases") + 
  geom_rect(aes(xmin=2000,ymin=0, xmax=2008, ymax=2630), fill = "black", alpha = 0.05) + 
  geom_ribbon(data = dat, aes(x = time, ymax=qtHig, ymin=qtLow, group = 1), fill="dodgerblue4", alpha =0.3) + 
  geom_line(data = dat, aes(x = time, y = median, group = 1, color="Training"), size = 0.9) + 
  geom_ribbon(data = datfut, aes(x = time, ymax=qtHig, ymin=qtLow, group = 1), fill="darkcyan", alpha =0.3) + 
  geom_line(data = datfut, aes(x = time, y = median, group = 1, color="Hindcast"), linetype = 1, size = 0.9) + 
  #geom_point(data = datfut, aes(x = time, y = median, group = 1), shape = 21, 
  #           colour = "darkcyan", fill = "darkcyan", size = 1, stroke = 1)  +
  geom_line(data = datL, aes(x = time, y = Pf, group = 1, color="Observed"), size = 0.9)  +
  geom_line(data= data.frame(x=2000,y=1:2650), aes(x=x,y=y), linetype = "dotted") + 
  geom_segment(aes(x = 2000.1, y = 2630, xend = 2007.95, yend = 2630), arrow = arrow(length = unit(0.35, "cm"))) + 
  geom_text(aes(x = 2001.3, 2480, label="Prediction"), size=4, fontface="bold") + 
  scale_color_manual(values = c(Training = "dodgerblue4", Hindcast = "darkcyan", Observed = "red2"))+
  scale_x_continuous(breaks = seq(from = 1975, to = 2007, by = 5), expand = c(0.02,0.02)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 2550, by = 500), expand = c(0.05,0.02), limits=c(0,2700)) +
  theme(strip.background = element_blank(), legend.justification = c(0,1), 
        legend.position = "none", legend.title= element_blank(), 
        legend.text = element_text(size = 9), 
        legend.background = element_rect(fill=NA),
        axis.text=element_text(size=12, color = "black"), 
        axis.title = element_text(size=13),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.line.x=element_line(),axis.line.y=element_line())

ggsave("figure2.pdf", p1, width = 12, height = 3.5)

