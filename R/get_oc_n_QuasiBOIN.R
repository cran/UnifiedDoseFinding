get_oc_QuasiBOIN <- function (target, p.true, score, ncohort, cohortsize, n.earlystop = 100,
                    startdose = 1, p.saf = 0.6 * target, p.tox = 1.4 * target,
                    cutoff.eli = 0.95,
                    extrasafe = FALSE, offset = 0.05, ntrial = 1000, seed = 100)
{
  if (target < 0.05) {
    stop("The target is too low! \n")
    return()
  }
  if (target > 0.6) {
    stop("The target is too high! \n")
    return()
  }
  if ((target - p.saf) < (0.1 * target)) {
    stop("The probability deemed safe cannot be higher than or too close to the target! \n")
    return()
  }
  if ((p.tox - target) < (0.1 * target)) {
    stop("The probability deemed toxic cannot be lower than or too close to the target! \n")
    return()
  }
  if (offset >= 0.5) {
    stop("The offset is too large! \n")
    return()
  }
  if (n.earlystop <= 6) {
    warning("The value of n.earlystop is too low to ensure good operating characteristics. Recommend n.earlystop = 9 to 18 \n")
    return()
  }
  set.seed(seed)
  ndose = nrow(p.true)
  npts = ncohort * cohortsize
  Y = matrix(rep(0, ndose * ntrial), ncol = ndose)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  dselect = rep(0, ntrial)
  temp = get_boundary_QuasiBOIN(target, ncohort, cohortsize, n.earlystop,
                      p.saf, p.tox, cutoff.eli, extrasafe, print = FALSE)
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  for (trial in 1:ntrial) {
    y <- rep(0, ndose)
    n <- rep(0, ndose)
    earlystop = 0
    d = startdose
    elimi = rep(0, ndose)
    for (i in 1:ncohort) {
      print(d)
      y[d] = y[d] + sum(as.numeric(t(rmultinom(cohortsize,1,prob =  p.true[d,]))%*%score))
      n[d] = n[d] + cohortsize
      if (n[d] >= n.earlystop)
        break
      if (!is.na(b.elim[n[d]])) {
        if (y[d] >= b.elim[n[d]]) {
          elimi[d:ndose] = 1
          if (d == 1) {
            earlystop = 1
            break
          }
        }
        if (extrasafe) {
          if (d == 1 && n[1] >= 3) {
            if (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                          1) > cutoff.eli - offset) {
              earlystop = 1
              break
            }
          }
        }
      }
      if (y[d] <= b.e[n[d]] && d != ndose) {
        if (elimi[d + 1] == 0)
          d = d + 1
      }
      else if (y[d] >= b.d[n[d]] && d != 1) {
        d = d - 1
      }
      else {
        d = d
      }
    }
    Y[trial, ] = y
    N[trial, ] = n
    if (earlystop == 1) {
      dselect[trial] = 99
    }
    else {
      dselect[trial] = select_mtd_QuasiBOIN(target, n, y, cutoff.eli,
                                  extrasafe, offset, print = FALSE)
    }
  }
  selpercent = rep(0, ndose)
  nptsdose = apply(N, 2, mean)
  ntoxdose = apply(Y, 2, mean)
  for (i in 1:ndose) {
    selpercent[i] = sum(dselect == i)/ntrial * 100
  }
  # cat("selection percentage at each dose level (%):\n")
  # cat(formatC(selpercent, digits = 1, format = "f"), sep = "  ",
  #     "\n")
  # cat("number of patients treated at each dose level:\n")
  # cat(formatC(nptsdose, digits = 1, format = "f"), sep = "  ",
  #     "\n")
  # cat("number of toxicity observed at each dose level:\n")
  # cat(formatC(ntoxdose, digits = 1, format = "f"), sep = "  ",
  #     "\n")
  # cat("average number of toxicities:", formatC(sum(Y)/ntrial,
  #                                              digits = 1, format = "f"), "\n")
  # cat("average number of patients:", formatC(sum(N)/ntrial,
  #                                            digits = 1, format = "f"), "\n")
  # cat("percentage of early stopping due to toxicity:", formatC(sum(dselect ==
  #                                                                    99)/ntrial * 100, digits = 1, format = "f"), "% \n")
  # if (length(which(p.true == target)) > 0) {
  #   if (which(p.true==target)==1){underdosing60 =underdosing80 =0}
  #   if (which(p.true==target)==2){underdosing60=mean(N[,1] > 0.6*npts)*100;underdosing80=mean(N[,1] > 0.8*npts)*100}
  #   if (which(p.true==target)>=3){
  #     if(dim(N)[1]>1){
  #       underdosing60=mean(rowSums(N[,1:(which(p.true==target)-1)])> 0.6*npts)* 100;underdosing80=mean(rowSums(N[,1:(which(p.true==target)-1)])> 0.8*npts)* 100;
  #     }else{
  #       underdosing60=mean(sum(N[,1:(which(p.true==target)-1)])> 0.6*npts)* 100;underdosing80=mean(sum(N[,1:(which(p.true==target)-1)])> 0.8*npts)* 100;
  #     }
  #   }
  #
  #
  #
  #   if (which(p.true == target) == ndose - 1) {
  #     overdosing60 = mean(N[, p.true > target] > 0.6 *
  #                           npts) * 100
  #     overdosing80 = mean(N[, p.true > target] > 0.8 *
  #                           npts) * 100
  #   }
  #   else {
  #     if(ntrial==1){
  #       overdosing60 = mean(sum(N[, p.true > target]) >
  #                             0.6 * npts) * 100
  #       overdosing80 = mean(sum(N[, p.true > target]) >
  #                             0.8 * npts) * 100
  #     }else{
  #       overdosing60 = mean(rowSums(N[, p.true > target]) >
  #                             0.6 * npts) * 100
  #       overdosing80 = mean(rowSums(N[, p.true > target]) >
  #                             0.8 * npts) * 100
  #     }
  #   }
  #   cat("risk of poor allocation:", formatC(mean(N[, p.true ==
  #                                                    target] < npts/ndose) * 100, digits = 1, format = "f"),
  #       "% \n")
  #   cat("risk of overdosing (>60% of patients treated above the MTD):",
  #       formatC(overdosing60, digits = 1, format = "f"),
  #       "% \n")
  #   cat("risk of overdosing (>80% of patients treated above the MTD):",
  #       formatC(overdosing80, digits = 1, format = "f"),
  #       "% \n")
  #   cat("risk of underdosing (>60% of patients treated under the MTD):",
  #       formatC(underdosing60, digits = 1, format = "f"),
  #       "% \n")
  #   cat("risk of underdosing (>80% of patients treated under the MTD):",
  #       formatC(underdosing80, digits = 1, format = "f"),
  #       "% \n")
  # }
  # if (length(which(p.true == target)) > 0) {
  #   list(target = target, p.true = p.true,
  #        ncohort = ncohort, cohortsize = cohortsize, startdose = startdose,
  #        p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli,
  #        extrasafe = extrasafe, offset = offset, ntrial = ntrial,
  #        dose = 1:ndose, selpercent = selpercent, nptsdose = nptsdose,
  #        ntoxdose = ntoxdose, totaltox = sum(Y)/ntrial, totaln = sum(N)/ntrial,
  #        pctearlystop = sum(dselect == 99)/ntrial * 100, overdose60 = overdosing60,
  #        overdose80 = overdosing80, underdose60 = underdosing60, underdose80=underdosing80)
  # }
  # else {
  #   list(target = target, p.true = p.true,
  #        ncohort = ncohort, cohortsize = cohortsize, startdose = startdose,
  #        p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli,
  #        extrasafe = extrasafe, offset = offset, ntrial = ntrial,
  #        dose = 1:ndose, selpercent = selpercent, nptsdose = nptsdose,
  #        ntoxdose = ntoxdose, totaltox = sum(Y)/ntrial, totaln = sum(N)/ntrial,
  #        pctearlystop = sum(dselect == 99)/ntrial * 100,overdose60 = overdosing60,
  #        overdose80 = overdosing80, underdose60 = underdosing60, underdose80=underdosing80)
  # }

  list(target = target, p.true = p.true,
       ncohort = ncohort, cohortsize = cohortsize, startdose = startdose,
       p.saf = p.saf, p.tox = p.tox,selpercent = selpercent, nptsdose = nptsdose,
       ntoxdose = ntoxdose)
}
