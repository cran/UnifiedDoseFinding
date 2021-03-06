# get.boundary <- function (target, ncohort, cohortsize, n.earlystop = 100, p.saf = 0.6 *
#                             target, p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
#                           offset = 0.05, print = TRUE)
# {
#   npts = ncohort * cohortsize
#   b.e = NULL
#   b.d = NULL
#   for (n in 1:npts) {
#     lambda1 = log((1 - p.saf)/(1 - target))/log(target *
#                                                   (1 - p.saf)/(p.saf * (1 - target)))
#     lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
#                                                            target)/(target * (1 - p.tox)))
#     cutoff1 = lambda1 * n
#     cutoff2 = lambda2 * n
#     b.e = c(b.e, cutoff1)
#     b.d = c(b.d, cutoff2)
#   }
#   boundaries = rbind(b.e, b.d)
#   return(boundaries)
# }
get_boundary_QuasiBOIN <- function (target, ncohort, cohortsize, n.earlystop = 100, p.saf = 0.6 *
                            target, p.tox = 1.4 * target, cutoff.eli = 0.95, extrasafe = FALSE,
                          offset = 0.05, print = TRUE)
{
  density1 <- function(p, n, m1, m2) {
    pbinom(m1, n, p) + 1 - pbinom(m2 - 1, n, p)
  }
  density2 <- function(p, n, m1) {
    1 - pbinom(m1, n, p)
  }
  density3 <- function(p, n, m2) {
    pbinom(m2 - 1, n, p)
  }
  # the following commented code delted by Mu
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
  npts = ncohort * cohortsize
  ntrt = NULL
  b.e = NULL
  b.d = NULL
  elim = NULL
  for (n in 1:npts) {
    lambda1 = log((1 - p.saf)/(1 - target))/log(target *
                                                  (1 - p.saf)/(p.saf * (1 - target)))
    lambda2 = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -
                                                           target)/(target * (1 - p.tox)))
    cutoff1 = lambda1 * n
    cutoff2 = lambda2 * n
    ntrt = c(ntrt, n)
    b.e = c(b.e, cutoff1)
    b.d = c(b.d, cutoff2)
    elimineed = 0
    if (n < 3) {
      elim = c(elim, NA)
    }
    else {
      for (ntox in 1:n) {
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) >
            cutoff.eli) {
          elimineed = 1
          break
        }
      }
      if (elimineed == 1) {
        elim = c(elim, ntox)
      }
      else {
        elim = c(elim, NA)
      }
    }
  }
  for (i in 1:length(b.d)) {
    if (!is.na(elim[i]) && (b.d[i] > elim[i]))
      b.d[i] = elim[i]
  }
  boundaries = rbind(ntrt, b.e, b.d, elim)[, 1:min(npts, n.earlystop)]
  rownames(boundaries) = c("Number of patients treated", "Escalate if # of DLT <=",
                           "Deescalate if # of DLT >=", "Eliminate if # of DLT >=")
  colnames(boundaries) = rep("", min(npts, n.earlystop))
  if (print) {
    message("Escalate dose if the observed toxicity rate at the current dose <= ",
        lambda1, "\n")
    message("Deescalate dose if the observed toxicity rate at the current dose >= ",
        lambda2, "\n\n")
    message("This is equivalent to the following decision boundaries\n")
    print(boundaries[, (1:floor(min(npts, n.earlystop)/cohortsize)) *
                       cohortsize])
    if (cohortsize > 1) {
      message("\n")
      message("A more completed version of the decision boundaries is given by\n")
      print(boundaries)
    }
    message("\n")
    if (!extrasafe)
      message("Default stopping rule: stop the trial if the lowest dose is eliminated.\n")
  }
  if (extrasafe) {
    stopbd = NULL
    ntrt = NULL
    for (n in 1:npts) {
      ntrt = c(ntrt, n)
      if (n < 3) {
        stopbd = c(stopbd, NA)
      }
      else {
        for (ntox in 1:n) {
          if (1 - pbeta(target, ntox + 1, n - ntox +
                        1) > cutoff.eli - offset) {
            stopneed = 1
            break
          }
        }
        if (stopneed == 1) {
          stopbd = c(stopbd, ntox)
        }
        else {
          stopbd = c(stopbd, NA)
        }
      }
    }
    stopboundary = rbind(ntrt, stopbd)[, 1:min(npts, n.earlystop)]
    rownames(stopboundary) = c("The number of patients treated at the lowest dose  ",
                               "Stop the trial if # of DLT >=        ")
    colnames(stopboundary) = rep("", min(npts, n.earlystop))
    if (print) {
      message("\n")
      message("In addition to the default stopping rule (i.e., stop the trial if the lowest dose is eliminated), \n")
      message("the following more strict stopping safety rule will be used for extra safety: \n")
      message(" stop the trial if (1) the number of patients treated at the lowest dose >= 3 AND",
          "\n", "(2) Pr(the toxicity rate of the lowest dose >",
          target, "| data) > ", cutoff.eli - offset, ",\n",
          "which corresponds to the following stopping boundaries:\n")
      print(stopboundary)
    }
  }
  if (!print)
    return(boundaries)
}
