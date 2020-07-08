get_oc_Quasi_CRM <- function (ptox, skeletons, target, score, cohortsize, ncohort, start.dose=1, mselection=1,
                              cutoff.eli	= 0.90, ntrial = 10, seed = 100)
{
  
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
  N <- matrix(rep(0, ntrial * dim(skeletons)[2]), ncol = dim(skeletons)[2])
  S <- matrix(rep(0, ntrial * dim(skeletons)[2]), ncol = dim(skeletons)[2])
  set.seed(seed)
  for(ii in 1:ntrial){
    ### define toxicity grade score system
    ngrade = length(score);
    s.max = max(score);
    target = target/s.max;  # standardize target ET score
    ndose = ncol(skeletons);
    y=rep(0, ndose);  #number of toxicity at each dose level
    n=rep(0, ndose);  #number of treated patients at each dose level
    dose.curr = start.dose;  # current dose level
    ptox.hat = numeric(ndose); # estimate of toxicity prob
    dose.select=rep(0, ndose); # a vector of indicators for dose selection
    stop=0; #indicate if trial stops early
    ### run a trial
    for(i in 1:ncohort)
    {
      # generate data for a new cohort of patients
      for(j in 1:cohortsize)
      {
        tox.ind = rmultinom1(1, 1, ptox[,dose.curr]);
        s.obs = sum(tox.ind*score)/s.max;
        y[dose.curr] = y[dose.curr] + s.obs;
      }

      n[dose.curr] = n[dose.curr] + cohortsize;

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

        # estimation based on the selected model
        p.overtox = integrate(posterior,lower=-Inf,upper=log(log(target)/log(skeletons[msel,1])), skeletons[msel,], y, n)$value/marginal[msel];
        if(p.overtox>cutoff.eli) { stop=1; break;}

        # calculate posterior mean of toxicity probability at each dose leavel
        for(j in 1:ndose) { ptox.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf, skeletons[msel,], y, n, j)$value/marginal[msel]; }
      }
      else   ### model averaging
      {
        pj.overtox = rep(0, nskel);
        for(k in 1:nskel)
        {
          pj.overtox[k] = integrate(posterior,lower=-Inf,upper=log(log(target)/log(skeletons[k,1])), skeletons[k,], y, n)$value/marginal[k];
        }
        p.overtox = sum(postprob*pj.overtox);
        if(p.overtox>cutoff.eli) { stop=1; break;}

        # calculate posterior mean of toxicity probability at each dose leavel
        ptoxj.hat = rep(0, nskel);
        for(j in 1:ndose){
          for(k in 1:nskel){
            ptoxj.hat[k] = integrate(posttoxf,lower=-Inf,upper=Inf, skeletons[k,], y, n, j)$value/marginal[k];
          }
          ptox.hat[j] = sum(postprob*ptoxj.hat);
        }
      }

      diff = abs(ptox.hat-target);
      dose.best = min(which(diff==min(diff)));
      #       dose.curr = dose.best;  # allowing dose skipping
      if(dose.best>dose.curr && dose.curr != ndose) dose.curr = dose.curr+1;
      if(dose.best<dose.curr && dose.curr != 1) dose.curr = dose.curr-1;
      
    }

    if(stop==0) { dose.select[dose.best] = dose.select[dose.best]+1; }
    N[ii, ] <- n
    S[ii, ] <- dose.select
    
  }
  selpercent <- apply(S, 2, mean) * 100
  nptsdose <- apply(N, 2, mean)
  list(selpercent = selpercent, nptsdose = nptsdose)
}