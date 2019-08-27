# Adapted code from MICE[1]
# [1] Buuren, S. V., & Groothuis-Oudshoorn, K. (2010). 
# mice: Multivariate imputation by chained equations in R. 
# Journal of statistical software, 1-68.
norm_draw <- function(y, ry, x, ridge = 1e-05, ...) {
    xobs <- x[ry, ]
    yobs <- y[ry]
    xtx <- crossprod(xobs)
    pen <- ridge * diag(xtx)
    if (length(pen) == 1) 
        pen <- matrix(pen)
    v <- solve(xtx + diag(pen))
    coef <- t(yobs %*% xobs %*% v)
    residuals <- yobs - xobs %*% coef
    df <- max(sum(ry) - ncol(x), 1)
    sigma.star <- sqrt(sum((residuals)^2)/rchisq(1, df))
    beta.star <- coef + (t(chol(v)) %*% rnorm(ncol(x))) * sigma.star
    parm <- list(coef, beta.star, sigma.star)
    names(parm) <- c("coef", "beta", "sigma")
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
  yobs <- y[ry]
  xtx <- t(xobs) %*% xobs

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