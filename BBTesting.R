# R script to test the BB solver in R.

require(BB)

resid <- function(p){
  r <- rep(NA, length(p))
  x <- p["x"]
  y <- p["y"]
#   r[1] <- x^2 - 8 - y
#   r[2] <- x - y - 2

  R1 <- x^2 - 8 - y
  R2 <- x - y - 2
#   
#   
#     with(as.list(p), {
#     r[1] <- x^2 - 8 - y
#     r[2] <- x - y - 2
#   })
# return(r)
  return(c(R1, R2))
}

p_init <- c(x=0, y=0)
ans <- BBsolve(par=p_init, fn=resid)