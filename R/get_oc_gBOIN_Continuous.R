get_oc_gBOIN_Continuous <-function(target,c_true,ncohort,cohortsize,ntrial,mu_1=0.6*target,mu_2=1.4*target,startdose = 1, seed = 100){
  set.seed(seed)
  ndose = length(c_true)
  npts = ncohort * cohortsize
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  for (trial in 1:ntrial){
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    d = startdose
    for (i in 1:ncohort) {
      y[d] = y[d] + sum(rnorm(cohortsize,c_true[d],0.05*d))#0.2*d))
      n[d] = n[d] + cohortsize
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
    dselect[trial] = select_mtd_gBOIN_Continuous(target, n, y)
  }
  selpercent = rep(0, ndose)
  nptsdose = apply(N, 2, mean)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  list(selpercent = selpercent, nptsdose = nptsdose)
}
