
## simulating dynamics
a <- 15 ## time it takes to consume W0 - W1 when canopy is closed
L0 <- 0.000001
T <- 40
l <- 1.5
C1 <- 0.8/365
C2 <- 0.2/365

calcG <- function(t) (t*(C1 + C2) - T*C2) / T

Lnext <- function(Lcurr){

  if(Lcurr >= 1) {
    Lnext <- 1
  }
  else {
    t <- a/Lcurr
    if (t > T) {
      G <- calcG(T)
    }
    else {
      G <- calcG(t)
    }
    Lnext <- Lcurr + 3.6*(G*T)^l
    if(Lnext > 1) Lnext <- 1
  }
    return(Lnext)
}


plot_L <- function(nT = 100) {

  data <- data.frame(T = 1:nT, L = rep(0, nT), dL = rep(0, nT), ddL = rep(0, nT))
  data[1,"L"] <- L0

  for (T in 2:nT) {

    data[T, "L"] <- Lnext(Lcurr = data[T-1, "L"])
    data[T, "dL"] <- data[T, "L"] - data[T-1, "L"]
    data[T, "ddL"] <- data[T, "dL"] - data[T-1, "dL"]

  }

  data <- data[data$T >1,]
  plot(data$T, data$L, pch = "")
  lines(data$T, data$L,)
  abline(v = 232)
  abline(v = 12)

  print(data)
  print(nrow(data[abs(data$ddL) > 1e-15, ]) / nrow(data))

}

plot_L(nT = 500)
