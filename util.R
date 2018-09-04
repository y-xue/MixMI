log_with_limits <- function(x,epsilon=1e-8) {
  # if (x < epsilon) {
  #   x = epsilon
  # }
  if (x == 0) {
  	x = .Machine$double.xmin
  }
  if (x > .Machine$double.xmax) {
    x = .Machine$double.xmax
  }
  return(log(x))
}