require(tgp)

#set.seed(1001)

s1 <- lhs(50000,rbind(
					c(0,0.5),  # sigOBS
					c(0,0.5),  # sigPRO
					c(0,100),  # MUEI
					c(0,100),  # MUIQ
					c(0,100),  # MUIS
					c(0,100),  # MUQS
					c(0,0.5),  # tau
					c(0,1),    # rho
					c(0,5),    # betaOUT
					 0.0,      # bT  
					c(0,10),   # bT4
					c(0,10),   # bT6 
					c(0,1),    # q0
					c(-10,10), # b1
					c(-10,10), # b2
					c(-10,10), # b3
					c(-10,10), # b4
					c(-10,10), # b5
					c(-10,10), # b6
					0.02     , # delta
					c(0.8,1),  # S_0
					c(0,0.1),  # E_0
					c(0,0.1),  # I_0
					c(0,0.1),  # Q_0
					c(0,1),    # K_0
					c(0,1)     # F_0
				)
			)



colnames(s1) <- c(
 					"sigOBS",
					"sigPRO",
					"muEI",		
					"muIQ",
					"muIS",
					"muQS",
					"tau",
					"rho",
					"betaOUT",
					"bT",
					"bT4",
					"bT6",
					"q0",
					"b1",
					"b2",
					"b3",
					"b4",
					"b5",
					"b6",
					"delta",
					"S_0",
					"E_0",
					"I_0",
					"Q_0",
					"K_0",
					"F_0"
				)

s1 <- as.data.frame(s1)
temp <- subset(s1,select=c("S_0","E_0","I_0","Q_0"))
temp2 <- t(apply(temp,1,function(x){x/sum(x)}))
s1[c("S_0","E_0","I_0","Q_0")] <- temp2

write.csv(round(s1,digits=4),file="parameters_TmeanB.csv",row.names=F)	# output rounded to 4 decimal digits

