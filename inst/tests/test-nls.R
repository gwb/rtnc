
require(rtnc)
require(numDeriv)


source('/Users/gwb/Hacks/Projects/BioAnalyzer/EXP_non_linear_regression/1-extract_data.R',chdir=T)


# rename variables
intensity <- obs
size <- s

mu0 <- sum(size*obs/sum(intensity))
sigma0 <- sqrt(sum(obs/sum(intensity)*(size-mu0)^2))
c0 <- max(intensity)/dnorm(0,0,sigma0) # kinda like a rescaling factor



fn <- function(x){
  
  f <- function(y){
      b0 <- y[1]
      b1 <- y[2]
      c <- y[3]
      u <- y[4]
      sigma <- y[5]
      res <- sum( (intensity - b0 - b1*size - c*dnorm(size, u, sigma))^2 )
      #if(is.na(res)){
      #  print(y)
      #}
      return(res)
  }

  b0 <- x[1]
  b1 <- x[2]
  c <- x[3]
  u <- x[4]
  sigma <- x[5]

  g.base = -2*(intensity - b0 - b1 * size - c*dnorm(size, u, sigma))
  g= c(0,0,0,0,0)

  g[1] = sum(g.base)
  g[2] = sum(g.base * size)
  g[3] = sum(g.base * dnorm(size, u, sigma))

  g[4] = sum((size-u)/sigma^2 * c*  dnorm(size, u, sigma) * g.base)
  g[5] = sum( g.base * (-c/sigma * dnorm(size, u, sigma) + (size-u)^2/sigma^3 * c * dnorm(size, u, sigma)))
  #og <- grad(f, x, method="simple")
  return(list(f(x), g))
}

nf <- function(x){
  f <- function(y) return((y[3]-3)^2 + 2)
  #g = c(0, 0, 2*(x[3] - 3))
  g = grad(f, x)
  return(list(f(x), g))
}

tnc(c(1,2,9000, mu0, sigma0), fn, low=c(NA,NA,0,100,0), up=c(10, NA, NA, NA, NA), maxnfeval=30, rescale=10)

tnc(c(-4,14,10), nf, low=c(NA, 0, 1.5))


est <- nls(intensity ~ b0 + b1*size + exp(logc)*dnorm(size, mu, exp(logsigma)),
           start=c(logc=log(c0), mu=mu0, logsigma=log(sigma0), b0=0, b1=0))

fn(c(0,0, c0, mu0, sigma0))
