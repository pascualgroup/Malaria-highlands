##################

simul <- "
	// compute transmission rate 
	double betaIN = exp(b1*season1+b2*season2+b3*season3+b4*season4+b5*season5+b6*season6+
	bT*tcovn2+bT4*Tm1t*season4+bT6*Tm2t*season6);

	// gamma white noise
	double dW = rgammawn(sigPRO,dt);

	// force of infection
	double foi = (betaOUT+((I+q0*Q)/pop)*betaIN)*dW/dt;
	
	// compute transition numbers
	double dBS = (delta*pop+dpopdt)*dt;
  double dSE = F*S*dt;
  double dEI = muEI*E*dt;
  double dIQ = muIQ*I*dt;
  double dIS = muIS*I*dt;
  double dQS = muQS*Q*dt;
  double dSD = delta*S*dt;
  double dED = delta*E*dt;
  double dID = delta*I*dt;
  double dQD = delta*Q*dt;
  double dKK = ((foi-K)/(tau/2.0))*dt;
  double dFF = ((K-F)/(tau/2.0))*dt;

	// compute equations
	S += dBS + dIS + dQS - dSE - dSD;
	E += dSE - dEI - dED;
	I += dEI - dIS - dIQ - dID;
	Q += dIQ - dQS - dQD;
  K += dKK;
  F += dFF;
	cases += rho*dEI;
	W += (dW-dt)/sigPRO;

	// track errors
	if (S < 0.0) { err -= S; S=0.0; }
	if (E < 0.0) { err -= E; E=0.0; }
	if (I < 0.0) { err -= I; I=0.0; }
	if (Q < 0.0) { err -= Q; Q=0.0; }
	if (K < 0.0) { err -= K; K=0.0; }
	if (F < 0.0) { err -= F; F=0.0; }
"

rmeas <- "
	double size = 1.0/sigOBS/sigOBS;
	double prob = size/(cases+size);
	Pf = rnbinom(size,prob);
"

dmeas <- "
	double size = 1.0/sigOBS/sigOBS;
	double prob = size/(cases+size);
	static double tol = 1.0e-18;
	lik = dnbinom(Pf,size,prob+tol,0)+tol;
	if (give_log) lik = log(lik);
"

fromEst <- Csnippet("
	TsigOBS = exp(sigOBS);
	TsigPRO = exp(sigPRO);
	TmuEI = exp(muEI);
	TmuIQ = exp(muIQ);
	TmuIS = exp(muIS);
	TmuQS = exp(muQS);
	Ttau = exp(tau);
	Trho = 1.0/(1.0+exp(-rho));
	TbetaOUT = exp(betaOUT);      
	Tq0 = 1.0/(1.0+exp(-q0));
	TbT4 = exp(bT4);
  TbT6 = exp(bT6);
	TS_0 = exp(S_0);
	TE_0 = exp(E_0);
	TI_0 = exp(I_0);
	TQ_0 = exp(Q_0);
	TK_0 = exp(K_0);
	TF_0 = exp(F_0);
	double sum = TS_0+TE_0+TI_0+TQ_0;
	TS_0 /= sum;
	TE_0 /= sum;
	TI_0 /= sum;
	TQ_0 /= sum;
")

toEst <- Csnippet("
	TsigOBS = log(sigOBS);
	TsigPRO = log(sigPRO);
	TmuEI = log(muEI);
	TmuIQ = log(muIQ);
	TmuIS = log(muIS);
	TmuQS = log(muQS);
	Ttau = log(tau);
	Trho = logit(rho);
	TbetaOUT = log(betaOUT);
	Tq0 = logit(q0);
	TbT4 = log(bT4);
	TbT6 = log(bT6);
	TK_0 = log(K_0);
	TF_0 = log(F_0);
	double sum = S_0+E_0+I_0+Q_0;
	TS_0 = log(S_0/sum);
	TE_0 = log(E_0/sum);
	TI_0 = log(I_0/sum);
	TQ_0 = log(Q_0/sum);
")

initlz <- Csnippet("
  double m = pop/(S_0+E_0+I_0+Q_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  Q = nearbyint(m*Q_0);
  K = K_0;
  F = F_0;
  cases = 0;
  W = 0;
  err = 0;
")

po <- pomp(
  data=subset(dat,select=c('time','Pf')),     
  times='time',
  t0=2*dat$time[1]-dat$time[2],
  covar=subset(dat,select=c('time','tcovn2','Tm1t','Tm2t','pop','dpopdt',sprintf('season%d',1:6))),
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

