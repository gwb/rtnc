\name{tnc}
\alias{tnc}
\title{Optimization based on the Truncated Newton Constrained algorithm}
\usage{
  tnc(x, fn, low = NULL, up = NULL, maxCGit = -1,
    maxnfeval = 20, eta = -1, stepmx = -1, accuracy = -1,
    fmin = 0, ftol = -1, xtol = -1, pgtol = -1,
    rescale = -1, coerce = T)
}
\arguments{
  \item{x}{: initial guess at argmin}

  \item{fn}{: a function that returns both value of the
  function to minimize and its gradient see details for
  more informations}

  \item{low,}{up : respectively lower and upper bounds of
  search space.  see details for more informations}

  \item{maxCGit}{: max. number of hessian * vector
  evaluation per main iteration if maxCGit == 0, the
  direction chosen is -gradient if maxCGit < 0, maxCGit is
  set to max(1, min(50, n/2))}

  \item{maxnfeval}{: max number of function evaluations}

  \item{eta}{: severity of the line search. if < 0 or > 1,
  set to 0.25}

  \item{stepmx}{: maximum step for the line search. May be
  increased during call.  if too small, will be set to
  10.0}

  \item{accuracy}{: relative precision for finite
  difference calculations if <= machine_precision, will be
  set to sqrt(machine_precision)}

  \item{fmin}{: minimum function value estimate}

  \item{ftol}{: precision goal for the value of f in the
  stopping criterion if ftol < 0.0, ftol is set to
  accuracy}

  \item{xtol}{: precision goal for the value of f in the
  stopping criterion if xtol < 0.0, xtol is set to
  sqrt(machine_precision)}

  \item{pgtol}{: precision goal for the value of the
  projected gradient in the stopping criterion (after
  applying x scaling factors) if pgtol < 0.0, pgtol is set
  to 1e-2 * sqrt(accuracy) setting it to 0.0 is not
  recommended}

  \item{rescale}{: f scaling factor (in log10) used to
  trigger f value rescaling if 0, rescale at each iteration
  if a big value, never rescale if < 0, rescale is set to
  1.3}

  \item{coerce}{: If true, force argument coercion to right
  type. Otherwise, throw an error if type is not correct.}
}
\value{
  x : argmin of function to minimize

  f : function evaluated at argmin

  g : gradient of function at argmin

  status : code returned by C underlying code see details
  for more informations

  nfeval : number of function evaluations
}
\description{
  Minimizes a function under constraints, using the TNC
  algorithm
}
\details{
  The function wraps a call to a C function written by
  Jean-Sebastien Roy. The R layer performs some simple type
  checks and coercions.

  The fn argument should be a function that takes a vector
  of length n and returns list of length two, whose first
  element is a scalar representing the value of the
  function to minimize evaluated at the argument vector,
  and the second element is a vector representing the
  gradient of the function evaluated at the argument. See
  examples.

  The arguments `low` and `up` indicate the minimum and
  maximum values of the parameter space. They should be
  vectors of size `length(x)`, where each index correspond
  to a bound on the corresponding dimension of the
  parameter space. If no constraint is desired on the this
  dimension, then low[dimension_i] should be set to NA. If
  no lower or upper constraint is desired, in any
  dimension, then the low and up vectors can be set to
  NULL. See examples.

  Possible values for 'status' are: -3 => Memory allocation
  failed (in C process) -2 => Invalid Parameters (n<0) -1
  => Infeasible (low bound > up bound) 0 => Local minima
  reach (|pg| ~= 0) 1 => Converged (|f_n - f_(n-1)| ~= 0) 2
  => Converged (|x_n - x_(n-1)| ~= 0) 3 => Max number of
  function evaluations reached 4 => Linear search failed 5
  => All lower bounds are equal to the upper bounds 6 =>
  User requested end of minimization
}
\examples{
# Defines the function to be minimized, and its gradient
fn <- function(x){
 x <- as.double(x)
 f = x[1]^2 +  abs(x[2])^2
 g = c(0,0)
 g[1] = 2 * x[1]
 g[2] = 3 * abs(x[2])^2
 g[2] = ifelse(x[2] < 0, -g[2], g[2])
 return(list(f, g))
}

# Finds the minimum of the function f(x) = x_1^2 + abs(x_2)^2, where x=(x_1,x_2)
# under the constraint that x_2 >= 1
tnc(c(1,2), fn, low=c(NA, 1))
}
\references{
  Stephen G. Nash (1984) "Newton Type Minimization Via the
  Lanczos Method", SIAM Journal of Numerical Analysis 21,
  pp. 770-778
}

