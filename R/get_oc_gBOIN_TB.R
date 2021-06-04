get_oc_gBOIN_TB<-function(target,pmat,weight,ncohort,cohortsize,n.earlystop = 100, ntrial,mu_1=0.6*target,mu_2=1.4*target,startdose = 1, seed = 100){
  if (n.earlystop <= 6) {
    cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n")
    return()
  }
  set.seed(seed)
  ndose = length(pmat)
  npts = ncohort * cohortsize
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  for (trial in 1:ntrial){
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    d = startdose
    score = 0
    for (i in 1:ncohort) {
      y0 =NULL
      for(irow in 1:5){
        y0=rbind(y0,rowSums(rmultinom(cohortsize,1,prob =  pmat[[d]][irow,])))
      }
      y[d] = y[d] + sum(y0*weight)
      n[d] = n[d] + cohortsize
      if (n[d] >= n.earlystop)
        break
      if (y[d]/n[d] <=(target+mu_1)/2&& d != ndose) {
          d = d + 1
      }
      else if (y[d]/n[d] >= (target+mu_2)/2 && d != 1) {
        d = d - 1
      }
      else {
        d = d
      }
    }
    Y[trial, ] = y
    N[trial, ] = n
    dselect[trial] = select_mtd_gBOIN_TB(target, n, y)
  }
  selpercent = rep(0, ndose)
  nptsdose = apply(N, 2, mean)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  list(selpercent = selpercent, nptsdose = nptsdose)
}
