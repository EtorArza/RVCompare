library(stats)


################## Comparison Functions ##################



#' The dominance rate of X_A over X_B given the density functions.
#'
#' Returns a real number in the interval [0,1] that represents the dominance rate of X_A over X_B.
#' Basically, we are measuring the amount of mass of X_A in which the cumulative distribution of X_A is higher minus the amount of mass of X_B in which the cumulative distribution of X_B is higher.
#  This value is then normalized so that all sections in which the distributions are equal are ignored, and finally, we apply the linear transformation $0.5(x -1)$ so that the dominance rate is defined between 0 and 1.
#'
#'
#' @param densityX_A The probability density function of the random variable X_A.
#' @param densityX_B The probability density function of the random variable X_B.
#' @param xlims an interval that represents the domain of definition the density functions.
#' @return Returns the dominance rate of X_A over X_B.
#' @seealso \code{\link{CpFromDensities}}
#' @export
#' @examples
### Example 1 ###
#' # If two symmetric distributions are centered in the same point (x = 0 in this case), then their Cd will be 0.5.
#' densityX_A <- normalDensity(0,1)
#' densityX_B <- uniformDensity(c(-2,2))
#' CdFromDensities(densityX_A, densityX_B, c(-5,5))
#'
#'
#' ### Example 2 ###
#' # If two distributions are equal, Cd will be 0.5.  Cd(X_A,X_A) = 0.5
#' CdFromDensities(densityX_A, densityX_A, c(-10,10))
#'
#'
#' ### Example 3 ###
#' # example on https://etorarza.github.io/pages/2021-interactive-comparing-RV.html
#' tau <- 0.11
#' densityX_A <- normalDensity(0.05,0.0015)
#' densityX_B <- mixtureDensity(c(normalDensity(0.05025,0.0015), normalDensity(0.04525, 0.0015)), xlims = c(.03,.07), weights = c(1 - tau, tau))
#' plot(densityX_A, from=0.03, to=0.07, type="l",  col="red", xlab="x", ylab="probability density")
#' curve(densityX_B, add=TRUE, col="blue", type="l", lty=2)
#' Cd <- CdFromDensities(densityX_A, densityX_B, c(.03,.07))
#' mtext(paste("Cd(X_A, X_B) =", format(round(Cd, 3), nsmall = 3)), side=3) # add Cd to plot as text
#' legend(x = c(0.0325, 0.045), y = c(200, 250),legend=c("X_A", "X_B"), col=c("red", "blue"), lty=1:2, cex=0.8) # add legend
#'
#'
#' ### Example 4 ###
#' densityX_A <- uniformDensity(c(0.1, 0.3))
#' densityX_B <- uniformDensity(c(-0.2,0.5))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-6,6))
#'
#' densityX_A <- mixtureDensity(c(uniformDensity(c(0.1,0.3)), uniformDensity(c(-1,-0.5))), xlims = c(-6,6))
#' densityX_B <- mixtureDensity(c(uniformDensity(c(-0.2,0.5)), uniformDensity(c(-1,-0.5))), xlims = c(-6,6))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-6,6))
CdFromDensities <- function(densityX_A, densityX_B, xlims, EPSILON = 1e-3) {

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
  }

  if(!isXlimsValid(xlims))
  {
    return(NULL)
  }


  if(!isFunctionDensity(densityX_A, xlims))
  {
    print("ERROR: argument densityX_A is not a non discrete probability density function.")
    return(NULL)
  }

  if(!isFunctionDensity(densityX_B, xlims))
  {
    print("ERROR: argument densityX_B is not a non discrete probability density function.")
    return(NULL)
  }

  cumX_A = cumulativeFromDensity(densityX_A, xlims)
  cumX_B = cumulativeFromDensity(densityX_B, xlims)
  cA = integrate(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_A(x)}, lower=xlims[[1]], upper=xlims[[2]])$value # the cA in the paper is cA^-1
  cB = integrate(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_B(x)}, lower=xlims[[1]], upper=xlims[[2]])$value # the cB in the paper is cB^-1

  if (min(cA, cB) < EPSILON) {
    return(0.5)
    print("WARNING: distributions are almost equal.")
  }

  f_to_integrate = function(y) { sapply(y, function(x) {
    if (cumX_A(x) > cumX_B(x) + EPSILON) {
      return(densityX_A(x) / cA)
    }else if(cumX_B(x) >  cumX_A(x) + EPSILON){
      return(-densityX_B(x) / cB)
    }else{
      return(0.0)
    }
    }) }

  return(0.5 *integrate(f_to_integrate, lower = xlims[[1]], upper = xlims[[2]])$value + 0.5)

}















#' The probability that X_A < X_B given the density functions.
#'
#' Returns a real number in the interval [0,1] that represents the probability that a sample observed from X_A is lower than a sample observed from X_B.
#'
#'
#'
#' @param densityX_A The probability density function of the random variable X_A.
#' @param densityX_B The probability density function of the random variable X_B.
#' @param xlims an interval that represents the domain of definition the density functions.
#' @return Returns the probability that X_A < X_B.
#' @seealso \code{\link{CdFromDensities}}
#' @export
#' @examples
#' ### Example 1 ###
#' # If two symmetric distributions are centered in the same point (x = 0 in this case), then their Cp will be 0.5.
#' densityX_A <- normalDensity(0,1)
#' densityX_B <- uniformDensity(c(-2,2))
#' Cp = CpFromDensities(densityX_A, densityX_B, c(-5,5))
#' plot(densityX_A, from=-5, to=5, type="l",  col="red", xlab="x", ylab="probability density")
#' curve(densityX_B, add=TRUE, col="blue", type="l", lty=2)
#' mtext(paste("Cp(X_A, X_B) =", format(round(Cp, 3), nsmall = 3)), side=3) # add Cp to plot as text
#' legend(x = c(-4.5, -2), y = c(0.325, 0.4),legend=c("X_A", "X_B"), col=c("red", "blue"), lty=1:2, cex=0.8) # add legend
#'
#'
#' ### Example 2 ###
#' # If two distributions are equal, Cp will be 0.5.  Cp(X_A,X_A) = 0.5
#' CpFromDensities(densityX_A, densityX_A, c(-10,10))
#'
#'
#' ### Example 3 ###
#' densityX_A <- normalDensity(-2,1)
#' densityX_B <- uniformDensity(c(-2,2))
#' # Cp(X_A,X_B) = 1 - Cp(X_B, X_A)
#' CpFromDensities(densityX_A, densityX_B, c(-Inf,Inf))
#' 1 - CpFromDensities(densityX_B, densityX_A, c(-Inf,Inf))
CpFromDensities <- function(densityX_A, densityX_B, xlims) {


  if(!isXlimsValid(xlims))
  {
    return(NULL)
  }

  if(!isFunctionDensity(densityX_A, xlims))
  {
    print("ERROR: argument densityX_A is not a non discrete probability density function.")
    return(NULL)
  }

  if(!isFunctionDensity(densityX_B, xlims))
  {
    print("ERROR: argument densityX_B is not a non discrete probability density function.")
    return(NULL)
  }

  f_to_integrate = function(y) { sapply(y, function(x) {densityX_A(x) * integrate(densityX_B, lower=x, upper=xlims[[2]])$value}) }

  return(integrate(f_to_integrate, lower = xlims[[1]], upper = xlims[[2]])$value)
}






################## Sanity Checks ##################
#' Check if a function is a (non-discrete) probability density function in a given domain.
#'
#' This function checks if an input function f is a non-discrete probability density function. For this to
#' be the case, the function needs to only return real values. The function also needs to be bounded, positive,
#' and its integral in the domain of definition needs to be 1.
#'
#'
#'
#' @param f the function to be checked.
#' @param xlims an interval that represents the domain of definition of f.
#' @param tol (optional parameter, default = 0.001) the integral of f is allowed to be in the interval (1-tol, 1+tol), to account for some reasonable error in the integration.
#' @return Returns True if the function is a non discrete probability density function. Otherwise, returns False.
#' @export
#' @examples
#' dist1 <- normalDensity(0,1)
#' isFunctionDensity(dist1, c(-2,2)) # the integral of the density of the normal distribution is too low in the interval (-2,2)
#' isFunctionDensity(dist1, c(-5,5)) # it is close enough from 1 in the interval (-5,5)

#' dist2 <- uniformDensity(c(0,1))
#' isFunctionDensity(dist2, xlims=c(-2,2))
#' isFunctionDensity(dist2, xlims=c(0.5,2)) # the integral is not 1
#'
#' dist3 <- function(x) 0.5/sqrt(x)
#' # The integral of the function being 1 is not enough to be considered a density function. It also needs to be boounded.
#' isFunctionDensity(dist3, c(1e-14,1))
#'
isFunctionDensity <- function(f, xlims, tol=1e-3) {

  if(tol > 0.1)
  {
    print(paste("ERROR: tol = ", toString(tol), " is not in the interval (0, 0.1)", sep=""))
    return(FALSE)
  }

  if (missing(f)) {
    print("ERROR: f not set.")
    return(FALSE)
  }


  if(!isXlimsValid(xlims))
  {
    return(FALSE)
  }

  # ignore bounded checks for infinite domain functions
  if (xlims[[1]] != -Inf && xlims[[2]] != Inf) {

    lowest_value = Inf
    arg_lowest_value = 0
    highest_value = -Inf
    arg_highest_value = 0

    nPointsChecked = 100
    for (i in 0:nPointsChecked) {
      x = xlims[[1]] + (xlims[[2]] - xlims[[1]]) * i / nPointsChecked
      if (f(x) < lowest_value) {
        lowest_value = f(x)
        arg_lowest_value = x
      }
      if(f(x) > highest_value)
      {
        highest_value = f(x)
        arg_highest_value = x
      }

    }



    if (lowest_value < 0 ) {
      print(paste("ERROR: for f to be correctly defined as a probability denstiy function, it needs to be positively defined in its domain. f(",toString(arg_lowest_value), ") = ", toString(lowest_value), sep = ""))
      return(FALSE)
    }
    if( highest_value > 1e6 / (xlims[[2]] - xlims[[1]]) )
    {
      print(paste("ERROR: the value of the density is too high in x = ", toString(arg_highest_value), ", where f(x) = ", toString(highest_value), ". This could mean that f is not bounded.", sep = ""))
      return(FALSE)
    }
  }

  integrand = integrate(f = f, lower = xlims[[1]], upper = xlims[[2]])["value"][[1]]
  if (  !(abs(integrand  -1 ) < tol)  )  {
    print(paste("ERROR: the integral of a density function in its domain must be 1. The value of the integral was ", toString(integrand), " instead.", sep=""))
    return(FALSE)
  }

  return(TRUE)
}








################## Probability density functions ##################

#' The probability density function of the normal distribution
#'
#' Returns the density function of the normal distribution with mean mu and standard deviation sigma.
#' The returned function is a single parameter function that returns the probability of the normal distribution in that point.
#' It is just a convinient wrapper of dnorm from the package 'stat' with some parameter checks.
#' @param mu the mean of the normal distribution.
#' @param sigma the standard deviation of the normal distribution.
#' @return Returns a callable function with a single parameter that describes the probability of the normal distribution in that point.
#' @family probability density distributions
#' @export
#' @examples
#' dist <- normalDensity(0,1)
#' dist(0)
normalDensity <- function(mu, sigma) {

  # parameter checks
  if (class(mu) != "numeric") {
    print("ERROR: mu is not numeric.")
    return(NULL)
  }
  if (class(sigma) != "numeric") {
    print("ERROR: sigma is not numeric.")
    return(NULL)
  }
  if (sigma <= 0) {
    print("ERROR: sigma needs to be positive.")
    return(NULL)
  }

  return(function(x) dnorm(x, mean = mu, sd = sigma))
}



#' The probability density function of the uniform distribution
#'
#' Returns the density function of the uniform distribution in the interval (xlims[[1]], xlims[[2]]).
#' The returned function is a single parameter function that returns the probability of the uniform distribution in that point.
#' It is just a convinient wrapper of dunif from the package 'stat' with some parameter checks.
#' @param xlims a tuple representing the interval of nonzero probability of the distribution.
#' @return Returns a callable function with a single parameter that rerturns the probability of the uniform distribution in each point.
#' @family probability density distributions
#' @export
#' @examples
#' dist <- uniformDensity(c(-2,2))
#' dist(-3)
#' dist(0)
#' dist(1)
uniformDensity <- function(xlims) {

  # parameter checks
  if (!isXlimsValid(xlims)) {
    print("ERROR: xlims in not correctly defined.")
    return(NULL)
  }

  return(function(x) dunif(x, min=xlims[[1]], max=xlims[[2]]))
}




#' A mixture of two or more distributions
#'
#' Returns the density function of the mixture distribution.
#' The returned function is a single parameter function that returns the probability of the mixture in that point.
#' @param densities the probability density functions to be combined.
#' @param xlims the domain of definition of the distributions in the mixture, and consequently, of the returned mixture density.
#' @param weights (optional) the weights of the distributions in the mixture. If it is not give, equal weights are assumed.
#' @return Returns a callable function with a single parameter that returns the probability of the mixture distribution each point.
#' @export
#' @examples
#' #If parameter weights not given, equal weights are assumed.
#' dist1 <- mixtureDensity(c(normalDensity(-2,1), normalDensity(2,1)))
#' plot(dist1, xlim = c(-5,5), xlab="x", ylab = "Probability density", main="Mixture of two Gaussians with equal weights", cex.main=0.85)
#'
#' dist2 <- mixtureDensity(c(normalDensity(-2,1), normalDensity(2,1)), weights=c(0.8,0.2))
#' plot(dist2, xlim = c(-5,5), xlab="x", ylab = "Probability density", main="Mixture of two Gaussians with different weights", cex.main=0.85)
mixtureDensity <- function(densities, xlims, weights=NULL) {


  # parameter checks
  if (length(densities) < 2) {
    print("ERROR: At least 2 distributions are required.")
    return(NULL)
  }

  n = length(densities)



  if (is.null(weights)) {
    weights = rep(1/n, n)
  }

  if (abs(sum(weights) - 1.0) > 1e-6) {
    print("ERROR: the sum of the weights must be 1.0.")
    return(NULL)
  }

  res <- function(x)
    {
      sum_of_probs = 0
      for (i in 1:n) {
        sum_of_probs = sum_of_probs + weights[[i]] * densities[[i]](x)
      }
      return(sum_of_probs)
    }
  return(res)
}

















################## Internal functions ##################

#' Check if xlims is a tuple that represents a valid interval in the real space.
#' @param xlims the tuple to be checked.
#' @return TRUE if it is a valid tuple. Otherwise prints error mesage and returns FALSE
#' @examples
#' xlims(c(-2,2))
#' xlims(c(2,-2))
isXlimsValid <- function(xlims) {

  if (missing(xlims)) {
    print("ERROR: xlims not set. Set the domain of definition with the parameter xlims. Example: xlims = c(-1,1) means the domain of definition is the interval (-1,1).")
    return(FALSE)
  }

  if (class(xlims) != "numeric") {
    print("ERROR: xlims is not numeric.")
    return(FALSE)
  }

  if (length(xlims) != 2) {
    print("ERROR: xlims is not a vector of length 2.")
    return(FALSE)
  }

  if (xlims[[1]] >= xlims[[2]]) {
    print("ERROR: xlims[[1]] needs to be lower than xlims [[2]] for xlims to be a valid interval.")
    return(FALSE)
  }

  return(TRUE)
}



#' Get the cumulative distribution function given the distribution function.
#' @param densityX The probability density function.
#' @param xlims the domain of definition of the density function.
#' @param sanityChecks (optional parameter, default = TRUE) boolean value indicating wether to check if the density function is correctly defined.
#' @return a callable function representing the cumulative distribution.
#' @keywords internal
cumulativeFromDensity <- function(densityX, xlims, sanityChecks = TRUE) {
  if(!isXlimsValid(xlims))
  {
    return(NULL)
  }
  if (sanityChecks) {
    if(!isFunctionDensity(densityX, xlims))
    {
      print("ERROR: argument densityX_A is not a non discrete probability density function.")
      return(NULL)
    }
  }



  return(  function(y) {sapply(y, function(x) {
    if(x < xlims[[1]] || x > xlims[[2]])
    {
      print(paste("ERROR: x = ", toString(x), " is out of the domain defined by xlims = ", toString(xlims), sep=""))
    }
    return(integrate(densityX, lower=xlims[[1]], upper=x)$value)
    }) }
  )
}







