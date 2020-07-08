next_mtd_BOIN <- function(target, n, y, d, p.saf = 0.6 * target, p.tox = 1.4 * target, cutoff.eli = 0.95,
                          extrasafe = FALSE, n.earlystop = 100){
  temp <- get_boundary_BOIN(target = target, ncohort = sum(n), cohortsize = 1, n.earlystop = n.earlystop,
                            p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli, extrasafe = extrasafe, print = FALSE)
  b.e = temp[2, ]
  b.d = temp[3, ]
  b.elim = temp[4, ]
  earlystop = 0
  ndose <- length(y)
  elimi = rep(0, ndose)
  if (!is.na(b.elim[n[d]])) {
    if (y[d] >= b.elim[n[d]]) {
      elimi[d:ndose] = 1
      if (d == 1) {
        earlystop = 1
        #break
      }
    }
    if (extrasafe) {
      if (d == 1 && n[1] >= 3) {
        if (1 - pbeta(target, y[1] + 1, n[1] - y[1] +
                      1) > cutoff.eli - offset) {
          earlystop = 1
          #break
        }
      }
    }
  }
  if (y[d] <= b.e[n[d]] && d != ndose) {
    if (elimi[d + 1] == 0&&n[d]>=2) d = d + 1 else d=d
  }
  else if (y[d] >= b.d[n[d]] && d != 1) {
    d = d - 1
  }
  else {
    d = d
  }
  if (earlystop == 1){
    99
  }
  else {
    d
  }

}


