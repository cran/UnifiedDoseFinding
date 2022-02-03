select_mtd_Ivanova_continuous <- function(target, c_resp, n){
  y <- sapply(c_resp, sum)
  le <- length(n)
  ## estimating the maximum tolerated dose
  trind=ifelse(n>0,1,0)*(1:le)
  if (n[1]==0) {trind=trind[-1]; dob=1} else dob=0
  #qstar=isot(c(y[trind],n[trind]))
  qstar=pava_continuous(y[trind]/n[trind],n[trind])
  ####  modify qstar so that the target dose is estimated according for the algorithm from Ivanova
  for (j in 1:length(qstar)) if (qstar[j]<target) qstar[j]=qstar[j]+1/1000/3^(7-j)
  list(dselect = order(round(abs(qstar-target),5))[1] + dob, n = n)
}
