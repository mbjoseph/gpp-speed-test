# generating data from full GP --------------------------------------------

# fixed values
n <- nmax
p <- 1
X <- matrix(1, nrow = n, ncol = 1)
coords <- data.frame(lon = runif(n), lat = runif(n))
D <- as.matrix(dist(coords))

# parameter values
eta <- .8
eta_sq <- eta ^ 2
sigma <- .3
phi <- 7
beta <- 2
C <- eta_sq * exp(-D * phi)
w <- c(t(chol(C)) %*% rnorm(n))

# response vector
log_mu <- c(X %*% beta + w + rnorm(n, sd = sigma))
y <- rpois(n, exp(log_mu))

# knot placement ----------------------------------------------------------

mx <- 4 # number of knots in x
my <- 4 # number of knots in y
m <- mx * my

place_knots <- function(mx, my, coords) {
  xcoords <- seq(min(coords[, 1]), max(coords[, 1]), length.out = mx + 1)
  ycoords <- seq(min(coords[, 2]), max(coords[, 2]), length.out = my + 1)
  x_offset <- diff(xcoords)[1] / 2
  y_offset <- diff(ycoords)[1] / 2
  expand.grid(lon = xcoords[-c(mx + 1)] + x_offset,
              lat = ycoords[-c(my + 1)] + y_offset)
}

coords_star <- place_knots(mx, my, coords)

D_star <- as.matrix(dist(coords_star))
D_site_star <- as.matrix(dist(rbind(coords, coords_star)))[1:n, (n + 1):(n + m)]


# bundle data for use in stan ---------------------------------------------

stan_d <- list(n = n,
               m = m,
               y = y,
               D_star = D_star,
               D_site_star = D_site_star,
               D = D,
               p = ncol(X),
               X = X)
