next_Ivanova_continuous <- function(target, eps, c_resp, n, d){
  y <- sapply(c_resp, sum)
  le <- length(n)
  qhat=y[d]/n[d]
  tstat=(qhat-target)/sd(c_resp[[d]])*sqrt(n[d])
  if (n[d]==1) d1=d
  if (n[d] >1) {	d1=d
  if (tstat< (-eps)) d1=d+1
  if (tstat> (eps)) d1=d-1
  }
  d=d1
  d=max(d,1)
  d=min(d,le)

  d

}
