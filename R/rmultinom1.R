rmultinom1 <- function(n, size, prob) {

  K <- length(prob) # #{classes}
  matrix(tabulate(sample(K, n * size, replace = TRUE, prob) + K * 0:(n - 1),
                  nbins = n * K),
         nrow = n, ncol = K, byrow = TRUE)
}
