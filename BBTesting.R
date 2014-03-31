# R script to test the BB solver in R.

require(BB)

resid <- function(p){
  with(as.list(p), {
    R1 <- x^2 - 8 - y
    R2 <- x - y - 2
    return(c(R1, R2))
  })
}

p_init <- c(x=0, y=0)
ans <- BBsolve(par=p_init, fn=resid)
print(ans$par)