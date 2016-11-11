#' Cissy Chan
#' Stochastic Simulation Project
#' M3S9
#' March 2014

set.seed(6)

# set up target function
f <- function(x) {
  (1/(x*(1-x))) * exp((-1/2)*(-2+log(x/(1-x)))^2)
}

# visualise the target function
x <- as.matrix(seq(0, 1, length=10000)[-c(1,10000)])
fx <- apply(x, 1, f)
plot(x, fx, type="l")

# fit the envelope function on f and find their intersections
# the envelope function part 1: exponential
e <- function(x) { exp(8 * 0.9 * (x-5/9)) }
ex <- apply(x, 1, e)
lines(x, ex, type="l")

# the envelope function part 2: flat line
lx <- rep(13.5, dim(x)[1])
lines(x, lx, type="l")

# the envelope function part 3: linear
h <- function(x) { -500 * x + 500 }
hx <- apply(x, 1, h)
lines(x, hx, type="l")

# intersection points of the three pieces of the envelope function
p <- log(13.5)/(7.2) + 5/9
q <- 486.5/500

# set up the full envelope function, i.e. "g(x)" unscaled
# scale this later when we calculate its area
g1 <- function(x) {
  ifelse(x <= p, exp(7.2 * (x-5/9)), 
         ifelse(x <= q, 13.5, -500 *x + 500))
}

# plot our unscaled envelope function
g1x <- apply(x, 1, g1)
plot(
  x=x, 
  y=fx, 
  type="l", 
  main="f(x) and g1(x)",
  xlab="x", 
  ylab="y"
)
lines(x, g1x, type="l", lwd=3)
legend("topleft", legend=c("g1(x)","f(x)"),lwd=c(3,1))

# area of g1 at intersections
s <- (1/7.2) * (exp(7.2 * (p-5/9))-exp(-4))
t <- s + 13.5 * (q-p)

# area under g1
area <- t + 250 * (q^2 - 2*q + 1)

# g; rescale g1 so that area = 1
g <- function(x) {
  g1(x)/area
}

# plot g(x) scaled
gx <- apply(x, 1, g)
lines(x, gx, type="l")

# check that g1 indeed envelopes f
optimise(function(k) f(k)-g1(k), c(0,1), maximum=T)

# using optimise() on the section of the function most likely to have the supremum
# get supremum M
M <- optimise(function(k) f(k)/g(k), c(0.5,1), maximum=T)$objective

# calculate theoretical acceptance probability
ap <- integrate(f, 0, 1)$value/M
ap

# set up CDF of g
G <- function(a) {
  ifelse(a <= p, (1/area) * (1/(8*0.9)) * (exp(8*0.9*a-4) - exp(-4)),
         ifelse(a <= q, (1/area) * (s + 13.5*(a-p)),
                (1/area) * (t + 250 * (q^2 - 2*q - a^2 + 2*a))))
}

Gx <- apply(x, 1, G)

# set up inverse of CDF of g
Ginv <- function(y) {
  ifelse(y <= G(p), (log(8*0.9*area*y + exp(-4))+4) / (8*0.9),
         ifelse(y <= G(q), (area*y-s) / 13.5 + p, 
                1 - sqrt(1 - ((area*y-t) / 250 - q^2 + 2*q))))
}

y <- x

Gy <- apply(y, 1, Ginv)

# plot G and its inverse
plot(x, Gx, type="l", main="G and G inverse", ylab="y")
lines(y, Gy, type="l", lwd=3)
legend("bottomright", c("G","G inverse"), lwd=c(1,3))

#---------------------------------------------------------------------------------------------------
# Rejection method
#---------------------------------------------------------------------------------------------------

# set up rejection algorithm
generate <- function(n) {
  rand <- vector('numeric', 0)
  m1 <- 0
  naccept <- 0
  total <- 0
  while (m1 < n) {
    count <- trunc((n-m1)*1.2)
    total <- total + count
    u <- runif(count)
    y <- Ginv(runif(count))
    naccept <- naccept + length( y[(M * g(y) * u) <= f(y)])
    rand <- c(rand, y[(M * g(y) * u) <= f(y)])
    m1 <- length(rand)
  }
  cat(paste0("number accepted = ", naccept, " out of ", total, "\n"))
  cat(paste0("acceptance prob = ", format(naccept/total),"\n"))
  rand[1:n]
}

# plot generated sample
hist(
  generate(1000000),
  breaks=150,
  freq=FALSE,
  main="Histogram for rejection algorithm",
  ylab="y",
  xlab="x"
)
lines(x, fx/2.506628)

#---------------------------------------------------------------------------------------------------
# Monte-Carlo integration
#---------------------------------------------------------------------------------------------------

# generate from uniform
monte1 <- function(n) {
  x <- runif(n)
  m <- (1/n) * sum(f(x)/1)
  var <- (1/n) * var(f(x)/1)
  result <- list(m=m, var=var)
}

# generate from envelope from before (Ginv function)
monte2 <- function(n) {
  r <- runif(n)
  x <- Ginv(r)
  m <- (1/n) * sum(f(x)/g(x))
  var <- (1/n) * var(f(x)/g(x))
  result <- list(m=m, var=var)
}

# generate samples from both methods
set.seed(6)
m1 <- monte1(1000000)
m2 <- monte2(1000000)

# variance reduction
m1$var
m2$var

# MC integration of f
m2$m

# stratified sampling: testing different functions for phi
phi1 <- function(x) { ifelse(x <= p,1,0) }
phi2 <- function(x) { ifelse(x > p & x <= q,1,0) }
phi3 <- function(x) { ifelse(x > q,1,0) }

ss1 <- function(n) {
  r <- runif(n)
  x <- Ginv(r)
  m <- (1/n) * sum(phi1(x) * f(x)/g(x))
  var <- (1/n) * var(phi1(x) * f(x)/g(x))
  result <- list(m=m, var=var)
}

ss2 <- function(n) {
  r <- runif(n)
  x <- Ginv(r)
  m <- (1/n) * sum(phi2(x) * f(x)/g(x))
  var <- (1/n) * var(phi2(x) * f(x)/g(x))
  result <- list(m=m, var=var)
}

ss3 <- function(n) {
  r <- runif(n)
  x <- Ginv(r)
  m <- (1/n)* sum(phi3(x) * f(x)/g(x))
  var <- (1/n) * var(phi3(x) * f(x)/g(x))
  result <- list(m=m, var=var)
}

# generate stratified samples
set.seed(5)
s1 <- ss1(100000) 
s2 <- ss2(100000) 
s3 <- ss3(100000)

s1$m + s2$m + s3$m
s1$var + s2$var + s3$var

#---------------------------------------------------------------------------------------------------
# Metropolis-Hastings sampler
#---------------------------------------------------------------------------------------------------

# set up algorithm
mh <- function(n_sim, burn_in, starting_point, v) {
  total <- 0
  # initialixe the chain
  X <- rep(starting_point, n_sim)
  for (i in 2:n_sim) {
    Y <- X[i-1] + rnorm(1,0,v)
    ifelse(Y>0 && Y<1, alpha <- f(Y)/f(X[i-1]), alpha <- 0)
    # generate uniform random sample for probability alpha occurring
    u <- runif(1)
    X[i] <- X[i-1] + (Y - X[i-1]) * (u<alpha)
    if (i>burn_in) {total <- total + as.numeric(u<alpha)}
  }
  accept_rate <- total / (n_sim - burn_in)
  X <- X[(burn_in+1):n_sim]
  result <- list(ar = accept_rate,x = X)
}

# how v affects acceptance rates
set.seed(1)
ar <- apply(as.matrix(seq(0.1,2.1,0.2)), 1, function(x) { mh(150,0,0.1,x)$ar })
ar2 <- apply(as.matrix(seq(0.1,2.1,0.2)), 1, function(x) { mh(150,0,0.1,x)$ar })
ar3 <- apply(as.matrix(seq(0.1,2.1,0.2)), 1, function(x) { mh(150,0,0.1,x)$ar })
ar4 <- apply(as.matrix(seq(0.1,2.1,0.2)), 1, function(x) { mh(150,0,0.1,x)$ar })
ar5 <- apply(as.matrix(seq(0.1,2.1,0.2)), 1, function(x) { mh(150,0,0.1,x)$ar })
cbind(v=seq(0.1, 2.1, 0.2), accept_rate=(ar+ar2+ar3+ar4+ar5)/5)

# plots showing convergence with different starting points
set.seed(5)
v <- 0.5
plot(
  seq(1, 150),
  mh(150, 0, 0.1, v)$x,
  type="l",
  main="MH algorithm convergence",
  xlab="n",
  ylab="X(n)"
)
lines(seq(1, 150), mh(150, 0, 0.5, v)$x, lwd=3)
lines(seq(1, 150), mh(150, 0, 0.9, v)$x, lwd=5)
legend("bottomright", c("x=0.1","x=0.5","x=0.9"), lwd=c(1,3,5))

# plots showing convergence with c = 0.3 to find burn-in number
set.seed(5)
c <- 0.3
plot(
  seq(1, 150),
  mh(150, 0, c, v)$x,
  type="l",
  main="MH algorithm convergence",
  xlab="n",
  ylab="X(n)"
)
lines(seq(1, 150), mh(150, 0, c, v)$x)
lines(seq(1, 150), mh(150, 0, c, v)$x)
lines(seq(1, 150), mh(150, 0, c, v)$x)
lines(seq(1, 150), mh(150, 0, c, v)$x)

# generate auto-correlation sequence
set.seed(5)
a <- mh(200000, 15, c, v)
a$ar
acf(a$x, main="Autocorrelation sequence for X, 1000 samples")
hist(a$x[1:200], breaks=50, freq=FALSE, main="MH samples 1", xlab="x", ylab="y")
lines(x, fx/2.5066)
hist(a$x[201:400], breaks=50, freq=FALSE, main="MH samples 2", xlab="x", ylab="y")
lines(x, fx/2.5066)
hist(a$x[401:600], breaks=50, freq=FALSE, main="MH samples 3", xlab="x", ylab="y")
lines(x, fx/2.5066)
hist(a$x[601:800], breaks=50, freq=FALSE, main="MH samples 4", xlab="x", ylab="y")
lines(x, fx/2.5066)
hist(a$x[801:1000], breaks=50, freq=FALSE, main="MH samples 5", xlab="x", ylab="y")
lines(x, fx/2.5066)
hist(a$x[1001:1200], breaks=50, freq=FALSE, main="MH samples 6", xlab="x", ylab="y")
lines(x, fx/2.5066)

# auto-correlation sequence after taking every 25 values
acf(
  a$x[seq(1, length(a$x), 25)],
  main="Autocorrelation sequence for x after taking every 25 values"
)

hist(
  a$x[seq(1, length(a$x), 25)],
  breaks=150, 
  freq=FALSE,
  main="MH random samples",
  xlab="x",
  ylab="y"
)

lines(x, fx/2.5066, lwd=3)


