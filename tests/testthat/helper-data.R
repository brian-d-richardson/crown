example_data <- function(seed = 42L, n = 1000L) {
  set.seed(seed)
  cluster <- rep(seq_len(20L), length.out = n)
  A <- rep(0:1, each = 10L)[cluster]
  S <- rep(c(1L, 0L), each = n / 2L)
  X <- rep(seq(-1, 1, length.out = 20L), length.out = n)
  W <- rnorm(n)
  R <- S * rbinom(n, 1, .8)
  C <- R * rbinom(n, 1, .2)
  Y <- rbinom(n, 1, plogis(-.3 + .2 * A + .4 * W + .1 * X))
  Y[R == 0 | C == 1] <- 0L
  data.frame(cluster, S, A, R, C, Y, X, W, wt = 1)
}
