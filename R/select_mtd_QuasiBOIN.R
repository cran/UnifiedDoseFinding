# select.mtd <- function (target, npts, ntox)
# {
#   pava <- function(x, wt = rep(1, length(x))) {
#     n <- length(x)
#     if (n <= 1)
#       return(x)
#     if (any(is.na(x)) || any(is.na(wt))) {
#       stop("Missing values in 'x' or 'wt' not allowed")
#     }
#     lvlsets <- (1:n)
#     repeat {
#       viol <- (as.vector(diff(x)) < 0)
#       if (!(any(viol)))
#         break
#       i <- min((1:(n - 1))[viol])
#       lvl1 <- lvlsets[i]
#       lvl2 <- lvlsets[i + 1]
#       ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
#       x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
#       lvlsets[ilvl] <- lvl1
#     }
#     x
#   }
#   y = ntox
#   n = npts
#   ndose = length(n)
#   adm.set = (n != 0)
#   adm.index = which(adm.set == T)
#   y.adm = y[adm.set]
#   n.adm = n[adm.set]
#   phat = (y.adm + 0.05)/(n.adm + 0.1)
#   phat = pava(phat, wt = n.adm)
#   phat = phat + (1:length(phat)) * 1e-10
#   selectd = sort(abs(phat - target), index.return = T)$ix[1]
#   selectdose = adm.index[selectd]
#   return(selectdose)
# }

select_mtd_QuasiBOIN <- function (target, npts, ntox, cutoff.eli = 0.95,
                        extrasafe = FALSE,
                        offset = 0.05, print = FALSE)
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
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= 3) {
      if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) >
          cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }
  if (extrasafe) {
    if (n[1] >= 3) {
      if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) >
          cutoff.eli - offset) {
        elimi[1:ndose] = 1
      }
    }
  }
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
    selectdose = 99
  }
  else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    phat = (y.adm + 0.05)/(n.adm + 0.1)
    phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm +
                                                           0.1)^2 * (n.adm + 0.1 + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10
    selectd = sort(abs(phat - target), index.return = T)$ix[1]
    selectdose = adm.index[selectd]
  }
  if (print == TRUE) {
    if (selectdose == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
    }
    else {
      message("The MTD is dose level ", selectdose, "\n\n")
    }
    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] -
                                 y[trtd] + 0.05))
    phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] +
                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                       (n[trtd] + 0.1 + 1))))
    message("Dose    Posterior DLT             95%                  \n",
        sep = "")
    message("Level     Estimate         Credible Interval   Pr(toxicity>",
        target, "|data)\n", sep = "")
    for (i in 1:ndose) {
      if (n[i] > 0) {
        message(" ", i, "        ", formatC(phat.all[i],
                                        digits = 2, format = "f"), "         (", formatC(qbeta(0.025,
                                                                                               y[i] + 0.05, n[i] - y[i] + 0.05), digits = 2,
                                                                                         format = "f"), ", ", formatC(qbeta(0.975, y[i] +
                                                                                                                              0.05, n[i] - y[i] + 0.05), digits = 2, format = "f"),
            ")            ", formatC(poverdose[i], digits = 2,
                                     format = "f"), "\n")
      }
      else {
        message(" ", i, "        ", "----", "         (",
            "------------", ")            ", "----", "\n")
      }
    }
    message("NOTE: no estimate is provided for the doses at which no patient was treated.")
  }
  else {
    return(selectdose)
  }
}
