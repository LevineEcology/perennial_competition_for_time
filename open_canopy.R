

calc_c <- function(t1, T, a_m = 0.3, gamma = 1, rr = 0.1) {

  return((1/T) * ((t1 * a_m) + (T-t1)*gamma*rr))

}

calc_G <- function(c, delta = 0.0815, nu = 1.5) {

  return((nu-1)/nu * c / delta)

}


mu_T <- 20
nint <- 100
T <- rexp(nint, rate = 1/mu_T)

B <- matrix(0, nrow = nint, ncol = nint)
N <- B
