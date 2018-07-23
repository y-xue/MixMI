norm_draw <- function(y, ry, x, ridge = 1e-05, ...) {
    xobs <- x[ry, ]
    yobs <- y[ry]
    xtx <- crossprod(xobs)  # SvB 21/01/2014
    pen <- ridge * diag(xtx)
    if (length(pen) == 1) 
        pen <- matrix(pen)
    v <- solve(xtx + diag(pen))
    coef <- t(yobs %*% xobs %*% v)
    residuals <- yobs - xobs %*% coef
    df <- max(sum(ry) - ncol(x), 1)  # SvB 31/10/2012
    sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))  # SvB 01/02/2011
    beta.star <- coef + (t(chol(v)) %*% rnorm(ncol(x))) * sigma.star
    parm <- list(coef, beta.star, sigma.star)  # SvB 10/2/2010
    names(parm) <- c("coef", "beta", "sigma")  # SvB 10/2/2010
    return(parm)
}

norm_fix <- function(y, ry, x, ridge=0.00001, ...)
{
# norm_fix
# Calculates regression coefficients + error estimate
#
# TNO Quality of Life
# authors: S. van Buuren and K. Groothuis-Oudshoorn
#
  xobs <- x[ry,]
  # print(xobs)
  yobs <- y[ry]
  xtx <- t(xobs) %*% xobs
  # print(dim(xobs))
  # print(dim(xtx))

  pen <- ridge * diag(xtx)
  if (length(pen)==1) pen <- matrix(pen)
  v <- solve(xtx + diag(pen))
  coef <- t(yobs %*% xobs %*% v)
  residuals <- yobs - xobs %*% coef
  sigma <- sqrt((sum(residuals^2))/sum(ry))
  parm <- list(coef, sigma)
  names(parm) <- c("beta", "sigma")
  return(parm)
}