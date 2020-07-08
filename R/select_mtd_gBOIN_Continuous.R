select_mtd_gBOIN_Continuous <- function (target, npts, ntox)
{
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
  y = ntox
  n = npts
  ndose = length(n)
  adm.set = (n != 0)
  adm.index = which(adm.set == T)
  y.adm = y[adm.set]
  n.adm = n[adm.set]
  phat = (y.adm + 0.05)/(n.adm + 0.1)
  phat = pava(phat, wt = n.adm)
  phat = phat + (1:length(phat)) * 1e-10
  selectd = sort(abs(phat - target), index.return = T)$ix[1]
  selectdose = adm.index[selectd]
  return(selectdose)
}
