get_oc_Ivanova_binary <- function(target, eps = 1, truetox,ncohort,cohortsize,n.earlystop = 100, ntrial,startdose=1, seed = 100){
  if (n.earlystop <= 6) {
    cat("Warning: the value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n")
    return()
  }
  set.seed(seed)
  ndose = length(truetox)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  for(trial in 1:ntrial){
    res = nonpart_binary(ga=target,truetox,ncohort,cohortsize,startdose, eps = eps, n.earlystop1 = n.earlystop)
    N[trial, ] = res$trials
    dselect[trial] = res$dselect
  }
  selpercent = rep(0, ndose)
  nptsdose = apply(N, 2, mean)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  list(selpercent = selpercent, nptsdose = nptsdose)
}
