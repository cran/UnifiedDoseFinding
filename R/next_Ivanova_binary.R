next_Ivanova_binary <- function(target, eps, y, n, d){
  le=length(y)
  qhat=y[d]/n[d]
  qhat1=max(qhat,0.001)
  tstat=(qhat-target)/sqrt(qhat1*(1-qhat1))*sqrt(n[d])
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
