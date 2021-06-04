next_gBOIN_TB <- function(target, n, y, d, mu_1 = 0.6 * target, mu_2 = 1.4 * target){
  ndose <- length(y)
  if (y[d]/n[d] <=(target+mu_1)/2&& d != ndose) {
    d = d + 1
  }
  else if (y[d]/n[d] >= (target+mu_2)/2 && d != 1) {
    d = d - 1
  }
  else {
    d = d
  }
  d
}
