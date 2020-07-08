nonpart_continuous=function(ga,truetox,ncohort,cohortsize,startdose=1, eps)
{
  # input data
  # model=xy[1]	# dose-response scenario
  n=ncohort*cohortsize		# total sample size
  #eps=1	# desgin parameter delta
  grsize=cohortsize	# cohort size
  # startup=0	# if start-up is desired, startup=1; otherwise no start-up
  # startupsize=round(log(.5)/log(1-ga)+.1,0)

  #
  # if (model==1)	{truetox=c(.05,.10,.20,.30,.50,.70)		#pr(tox)
  # }
  # if (model==2)	{truetox=c(.30,.4,.52,.61,.76,.87)		#pr(tox)
  # }
  # if (model==3)	{truetox=c(.05,.06,.08,.11,.19,.34)		#pr(tox)
  # }
  # if (model==4)	{truetox=c(.06,.08,.12,.18,.40,.71)		#pr(tox)
  # }
  # if (model==5)	{truetox=c(.00,.0,.03,.05,.11,.22)		#pr(tox)
  # }

  le=length(truetox)
  trials=rep(0,le)
  y=rep(0,le)   		# toxicity outcome, y=1 if toxicity, y=0 if not
  dose0=startdose			#starting dose
  dosenum=dose0
  count=0
  c_true = truetox
  # record response for each dose under each trial
  c_resp = list()
  for(i in 1:le) c_resp[[i]] = 0

  # if (startup==1) {while (sum(y)==0 && sum(trials)<startupsize*(le-dose0+1) ) {
  # trials[dosenum]=trials[dosenum]+startupsize
  # resp=rbinom(1,startupsize,truetox[dosenum])
  # y[dosenum]=y[dosenum]+resp
  # dosenum=ifelse(resp>0,max(dosenum-1,1),min(le,dosenum+1))
  # }
  # }

  count=sum(trials)
  while (count<n)	{
    groupsize=min(n-count,grsize)
    trials[dosenum]=trials[dosenum]+groupsize
    #resp=rbinom(1,groupsize,truetox[dosenum])
    resp = rnorm(cohortsize,c_true[dosenum],c_true[dosenum])#c_true[dosenum])#0.1*dosenum)
    c_resp[[dosenum]] = c(c_resp[[dosenum]],resp)
    y[dosenum]=y[dosenum]+sum(resp)
    qhat=y[dosenum]/trials[dosenum]
    tstat=(qhat-ga)/sd(c_resp[[dosenum]][-1])*sqrt(trials[dosenum])
    if (trials[dosenum]==1) dosenum1=dosenum
    if (trials[dosenum] >1) {	dosenum1=dosenum
    if (tstat< (-eps)) dosenum1=dosenum+1
    if (tstat> (eps)) dosenum1=dosenum-1
    }
    dosenum=dosenum1
    dosenum=max(dosenum,1)
    dosenum=min(dosenum,le)
    count=count+groupsize
  }

  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  ## estimating the maximum tolerated dose
  trind=ifelse(trials>0,1,0)*(1:le)
  if (trials[1]==0) {trind=trind[-1]; dob=1} else dob=0
  #qstar=isot(c(y[trind],trials[trind]))
  qstar=pava(y[trind]/trials[trind],trials[trind])
  ####  modify qstar so that the target dose is estimated according for the algorithm from Ivanova
  for (j in 1:length(qstar)) if (qstar[j]<ga) qstar[j]=qstar[j]+1/1000/3^(7-j)
  list(dselect = order(round(abs(qstar-ga),5))[1] + dob, trials = trials)
}
