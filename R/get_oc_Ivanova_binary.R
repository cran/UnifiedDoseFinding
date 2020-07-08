get_oc_Ivanova_binary <- function(target, eps = 1, truetox,ncohort,cohortsize,ntrial,startdose=1, seed = 100){
  set.seed(seed)
  ndose = length(truetox)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  for(trial in 1:ntrial){
    res = nonpart_binary(ga=target,truetox,ncohort,cohortsize,startdose, eps = eps)
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
