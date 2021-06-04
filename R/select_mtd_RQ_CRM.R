select_mtd_RQ_CRM <- function(target, n, y, score, skeleton, mselection = 1){
  # if a single skeleton is inputed as a vector, convert it to a matrix
  if(is.vector(skeletons)) skeletons=t(as.matrix(skeletons));

  nskel = nrow(skeletons);
  mprior = rep(1/nskel, nskel);  # prior for each model formed the skeleton


  # posterior = likelihood x prior
  posterior <- function(alpha, p, y, n)
  {
    sigma2 = 2;
    lik=1;
    for(j in 1:length(p))
    {
      pj = p[j]^(exp(alpha));
      lik = lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
  }


  # the posterior mean of ptox
  posttoxf <- function(alpha, p, y, n, j) { p[j]^(exp(alpha))*posterior(alpha, p, y, n); }


  s.max = max(score);
  target1 = target/s.max;  # standardize target ET score
  ndose = ncol(skeletons);
  ptox.hat = numeric(ndose); # estimate of toxicity prob
  dose.select=rep(0, ndose); # a vector of indicators for dose selection

  marginal = rep(0, nskel);
  for(k in 1:nskel)
  {
    marginal[k] = integrate(posterior,lower=-Inf,upper=Inf, skeletons[k,], y, n)$value;
  }

  postprob = (marginal*mprior)/sum(marginal*mprior);


  if(mselection==1)  ### model selection
  {
    # model selection, identify the model with the highest posterior prob
    if(nskel>1) { msel = which(postprob==max(postprob)); }
    else msel = 1;
    # calculate posterior mean of toxicity probability at each dose leavel
    for(j in 1:ndose) { ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, skeletons[msel,], y, n, j)$value/marginal[msel]; }
  }

  else   ### model averaging
  {
    # calculate posterior mean of toxicity probability at each dose leavel
    ptoxj.hat = rep(0, nskel);
    for(j in 1:ndose){
      for(k in 1:nskel){
        ptoxj.hat[k] = integrate(posttoxf,lower=-Inf,upper=Inf, skeletons[k,], y, n, j)$value/marginal[k];
      }
      ptox.hat[j] = sum(postprob*ptoxj.hat);
    }
  }

  diff = abs(ptox.hat-target1);
  dose.best = min(which(diff==min(diff)));

  dose.select[dose.best] = dose.select[dose.best]+1;


  which(dose.select == 1)



}
