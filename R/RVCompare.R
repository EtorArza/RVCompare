library(stats)
library(pracma)


################## Main functions ##################


#' Estimate X'_A and X'_B bounds with bootstrap 0.632
#'
#' Estimate the confidence intervals for the density estimations of X'_A and X'_B using bootstrap.
#' As a bonus, in addition to the density, the bounds of the cumulative density are also compared.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B.
#' @param nOfQuantiles the number of points in the interval [0,1] in which the density is estimated.
#' @param nOfBootstrapSamples (optional, default value 1e6) how many bootsrap samples to average.
#' @param alpha (optional, default value 0.2) the error of the confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param returnDataframe (optional, default value FALSE, as it returns a list by default) Wether to return a dataframe or a list.
#' @return Returns a list with the following fields:
#'
#' - p: values in the interval [0,1] that represent the nOfQuantiles points in which the densities are estimated. Useful for plotting.
#'
#' - X_prima_A_density_estimation: an array with the estimated probability densites of X_prima_A for each point (p[[i]] + p[[i+1]])/2.
#'
#' - X_prima_A_density_upper: an array with the upper bounds of confidence 1 - alpha of the density estimation of X_prima_A
#'
#' - X_prima_A_density_lower: an array with the lower bounds of confidence 1 - alpha of the density estimation of X_prima_A
#'
#' - X_prima_B_density_estimation: The same as X_prima_A_density_estimation for X'_B.
#'
#' - X_prima_B_density_upper: The same as X_prima_A_density_upper for X'_B
#'
#' - X_prima_B_density_lower: The same as X_prima_A_density_lower for X'_B
#'
#' - X_prima_A_cumulative_estimation, X_prima_A_cumulative_lower, X_prima_A_cumulative_upper, X_prima_B_cumulative_estimation, X_prima_B_cumulative_lower, X_prima_B_cumulative_upper: the same as the above fields, but it describes the cumulative distribution instead of the density.
#' @export
#' @examples
### Example 1 ###
#' X_A_observed <- c(0,1,1,2,2,3,4)
#' X_B_observed <- c(0,1,1,2,3,1,1,5,6)
#' res <- get_X_prima_AB_bounds_bootstrap(X_A_observed, X_B_observed, 100, returnDataframe=TRUE, nOfBootstrapSamples=1e3)
#'
#'
#' densityesPlot = ggplot() +
#'  geom_line(data = res, aes(x=p, y=X_prima_A_cumulative_estimation, colour = "X'_A", linetype="X'_A")) +
#'  geom_line(data = res, aes(x=p, y=X_prima_B_cumulative_estimation, colour = "X'_B",  linetype ="X'_B")) +
#'  scale_colour_manual("", breaks = c("X'_A", "X'_B"),  values = c("#F8766D", "#00BFC4")) +
#'  scale_linetype_manual("", breaks = c("X'_A", "X'_B"), values = c("dashed", "solid")) +
#'  xlab('x') +
#'  ylab('cumulative probability') +
#'  theme_bw()  # Black and white theme
#'  print(densityesPlot)
get_X_prima_AB_bounds_bootstrap <- function(X_A_observed, X_B_observed, nOfQuantiles, nOfBootstrapSamples=1e4, alpha=0.2,  EPSILON=1e-20, returnDataframe=FALSE) {

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
  }

  n <- length(X_A_observed)
  m <- length(X_B_observed)

  nDatapointsWhereDensityEstimated <- nOfQuantiles - 1

  p <- 0:(nDatapointsWhereDensityEstimated-1) / (nDatapointsWhereDensityEstimated - 1) # the size of the intervals is 1 / (nDatapointsWhereDensityEstimated - 1).

  dataA <- matrix(0, nrow = nOfBootstrapSamples, ncol = nDatapointsWhereDensityEstimated)
  dataB <- matrix(0, nrow = nOfBootstrapSamples, ncol = nDatapointsWhereDensityEstimated)

  pb = txtProgressBar(min = 1, max = nOfBootstrapSamples, initial = 1, style = 3)

  for (i in 1:nOfBootstrapSamples) {
    setTxtProgressBar(pb,i)
    bootStrapSampleA <- sample(X_A_observed, size= min(n,m), replace=TRUE)
    bootStrapSampleB <- sample(X_B_observed, size=min(n,m), replace=TRUE)

    ranksObj <- ranksOfObserved(bootStrapSampleA, bootStrapSampleB, EPSILON)

    X_A_ranks <- sort(ranksObj$X_A_ranks)
    X_B_ranks <- sort(ranksObj$X_B_ranks)
    r_max <- ranksObj$r_max
    j_max <- nDatapointsWhereDensityEstimated-1

    for (j in 0:j_max) {
      dataA[[i,j+1]] <- helperGet_X_prima_AB_bounds_bootstrap(sortedRanks = X_A_ranks, r_max = r_max, j = j, j_max = j_max)
      dataB[[i,j+1]] <- helperGet_X_prima_AB_bounds_bootstrap(sortedRanks = X_B_ranks, r_max = r_max, j = j, j_max = j_max)
    }
  }
  res <- list()

  res$p <- p

  quantiles <- apply(dataA, 2, quantile, probs = c(alpha/2, 0.5, 1.0 - alpha/2))
  res$X_prima_A_density_estimation <- quantiles[2,]
  res$X_prima_A_density_upper <- quantiles[1,]
  res$X_prima_A_density_lower<- quantiles[3,]

  res$X_prima_A_cumulative_estimation = cumsum(c(0,head(res$X_prima_A_density_estimation,-1)) / j_max)
  res$X_prima_A_cumulative_upper = cumsum(c(0,head(res$X_prima_A_density_upper,-1)) / j_max)
  res$X_prima_A_cumulative_lower = cumsum(c(0,head(res$X_prima_A_density_lower,-1)) / j_max)




  quantiles <- apply(dataB, 2, quantile, probs = c(alpha/2, 0.5, 1.0 - alpha/2))
  res$X_prima_B_density_estimation <- quantiles[2,]
  res$X_prima_B_density_upper <- quantiles[1,]
  res$X_prima_B_density_lower<- quantiles[3,]

  res$X_prima_B_cumulative_estimation = cumsum(c(0,head(res$X_prima_B_density_estimation,-1)) / j_max)
  res$X_prima_B_cumulative_upper = cumsum(c(0,head(res$X_prima_B_density_upper,-1)) / j_max)
  res$X_prima_B_cumulative_lower = cumsum(c(0,head(res$X_prima_B_density_lower,-1)) / j_max)

  if (returnDataframe) {
    df <- data.frame(matrix(unlist(res), nrow=length(res$p), byrow=FALSE))
    colnames(df) <- names(res)
    return(df)
  }else{
    return(res)
  }

}



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
#' densityX_B <- mixtureDensity(c(normalDensity(0.05025,0.0015), normalDensity(0.04525, 0.0015)), weights = c(1 - tau, tau))
#' plot(densityX_A, from=0.03, to=0.07, type="l",  col="red", xlab="x", ylab="probability density")
#' curve(densityX_B, add=TRUE, col="blue", type="l", lty=2)
#' Cd <- CdFromDensities(densityX_A, densityX_B, c(.03,.07))
#' mtext(paste("Cd(X_A, X_B) =", format(round(Cd, 3), nsmall = 3)), side=3) # add Cd to plot as text
#' legend(x = c(0.0325, 0.045), y = c(200, 250),legend=c("X_A", "X_B"), col=c("red", "blue"), lty=1:2, cex=0.8) # add legend
#'
#'
#' ### Example 4 ###
#' # The dominance factor ignores the mass of the probability where the distribution functinos are equal.
#' densityX_A <- uniformDensity(c(0.1, 0.3))
#' densityX_B <- uniformDensity(c(-0.2,0.5))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-2,2))
#'
#' densityX_A <- mixtureDensity(c(uniformDensity(c(0.1,0.3)), uniformDensity(c(-1,-0.5))))
#' densityX_B <- mixtureDensity(c(uniformDensity(c(-0.2,0.5)), uniformDensity(c(-1,-0.5))))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-2,2))
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
  cA = integral(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_A(x)}, xmin=xlims[[1]], xmax=xlims[[2]], method = "Simpson") # the cA in the paper is cA^-1
  cB = integral(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_B(x)}, xmin=xlims[[1]], xmax=xlims[[2]], method = "Simpson") # the cB in the paper is cB^-1

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

  return(0.5 *integral(f_to_integrate, xmin = xlims[[1]], xmax = xlims[[2]], method = "Simpson") + 0.5)

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
#' CpFromDensities(densityX_A, densityX_B, c(-8,4))
#' 1 - CpFromDensities(densityX_B, densityX_A, c(-8,4))
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


  integrand = integral(f, xmin = xlims[[1]], xmax = xlims[[2]])
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
mixtureDensity <- function(densities, weights=NULL) {


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

#' Check if xlims is a tuple that represents a valid bounded interval in the real space.
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

  if (xlims[[1]] == -Inf | xlims[[2]] == Inf) {
    print("ERROR: xlims needs to be a bounded interval.")
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
    return(integral(densityX, xmin=xlims[[1]], xmax=x, method="Simpson"))
    }) }
  )
}

#' Get the ranks from the values of observed X_A and X_B. Ranks go from 0 to r_max, where r_max is the number of unique values in c(X_A_observed, X_B_observed)
#' @param X_A_observed array of the samples (real values) of X_A.
#' @param X_B_observed array of the samples (real values) of X_B.
#' @param EPSILON (optional, default value 1e-20) when will two values be different.
#' @return a list with three fields: X_A_ranks, X_B_ranks and r_max (the number of unique values minus 1).
#' @keywords internal
ranksOfObserved <- function(X_A_observed, X_B_observed, EPSILON=1e-20) {

  n = length(X_A_observed)
  m = length(X_B_observed)

  all_values = c(X_A_observed, X_B_observed)

  order = order(all_values)
  order_original <- order
  inv_order = array(data=0, dim = length(all_values))

  for (i in 1:length(order)) {
    inv_order[order[[i]]] = i
  }


  n_repeated <- 0
  for (i in 1:(length(order_original)-1)) {
    if (abs(all_values[[order[[i]]]] - all_values[[order[[i+1]]]]) < EPSILON) {
      inv_order[[order[[i+1]]]] = inv_order[[order[[i]]]]
      n_repeated = n_repeated + 1
    }else{
      inv_order[[order[[i+1]]]] = inv_order[[order[[i+1]]]] - n_repeated
    }
  }



  # print(all_values)
  # print("order_original")
  # print(order_original)
  # print("order")
  # print(order)
  # print("inv_order")
  # print(inv_order)

  return(list("X_A_ranks"=inv_order[1:n]-1, "X_B_ranks"=inv_order[n+1:m]-1, "r_max"= n+m - n_repeated - 1))
}



#' Get a mixture of uniform (AKA tophat) distributions of size 'kernelSize' centered in the points 'kernelPositions'.
#' @param kernelPositions an array of real values that decribe the positions of the kernels.
#' @param kernelSize The size of the uniform kernel: a positive float.
#' @return a single parameter callable function that computes the probability density of the mixture in any point.
#' @keywords internal
mixtureOfUniforms <- function(kernelPositions, kernelSize) {

  n = length(kernelPositions)
  distributions = c()

  for (i in 1:n) {
    distributions = c(distributions,uniformDensity(c(kernelPositions[[i]] - kernelSize/2, kernelPositions[[i]] + kernelSize/2)))
  }

  return(mixtureDensity(distributions))
}



#' Get X'_A and X'_B from the observed values of X_A and X_B.
#'
#' X'_A and X'_B are two random variables defined in the interval [0,1] that have
#' the same C_p and C_d as the kernel density estimates of the observed X_A and X_B,
#' if a sufficiently small uniform kernel is used.
#' @param X_A_observed array of the samples (real values) of X_A.
#' @param X_B_observed array of the samples (real values) of X_B.
#' @param EPSILON (optional, default value 1e-20) when will two values be different.
#' @return a list with two fields X_prima_A and X_prima_B: each representing the density of the estimated distributions.
#' @export
#' @keywords internal
#' @examples
#' X_A_observed <- c(0,2,1)
#' X_B_observed <- c(1,6,1,3)
#' res <- get_X_prima_AB_density(X_A_observed, X_B_observed)
#' x = 0:1001/1001
#' matplot(x,cbind(res$X_prima_A(x),res$X_prima_B(x)),type="l",col=c("red","blue"), ylab='Probability density')
#' legend(x = c(0.7, 1.0), y = c(2.0, 2.5),legend=c("X'_A", "X'_B"), col=c("red", "blue"), lty=1:2, cex=0.8) # add legend
get_X_prima_AB_density <- function(X_A_observed, X_B_observed, EPSILON=1e-20) {

  ranksObj <- ranksOfObserved(X_A_observed, X_B_observed, EPSILON)

  X_A_ranks <- ranksObj$X_A_ranks
  X_B_ranks <- ranksObj$X_B_ranks
  r_max <- ranksObj$r_max

  kernelPositionsA <- X_A_ranks / (r_max + 1) + (0.5 / (r_max+1))
  kernelPositionsB <- X_B_ranks / (r_max + 1) + (0.5 / (r_max+1))

  kernelSize <- 1 / (r_max + 1)

  return(list("X_prima_A"=mixtureOfUniforms(kernelPositionsA, kernelSize), "X_prima_B"=mixtureOfUniforms(kernelPositionsB, kernelSize)))
}


#' Helper function for get_X_prima_AB_bounds_bootstrap.
#'
#' The density corresponding to the position in index j is computed, given the SORTED ranks, and r_max
#' @param sortedRanks the sorted ranks of either the observed X_A or X_B.
#' @param r_max The largest rank.
#' @param j the index that corresponds to the position. j goes from 0 to nDatapointsWhereDensityEstimated-1
#' @param j_max the largest index that will be used. its value is nDatapointsWhereDensityEstimated-1
#' @keywords internal
#' @examples
#' j_max <- 1000
#' f <- function(x){helperGet_X_prima_AB_bounds_bootstrap(sortedRanks=c(0,0,0,0,0,0,2,3,4), r_max=4, j=x, j_max=j_max)}
#' plot(x = 0:j_max / j_max, y = sapply(0:j_max, f), type="l")
#' plot(x = 0:j_max / j_max, y = cumsum(c(0,head(sapply(0:j_max, f),-1)) / j_max), type="l")
#' @export
#' @return the probability density in this point
helperGet_X_prima_AB_bounds_bootstrap <- function(sortedRanks, r_max, j, j_max) {
    p_corresponding_to_j <- j / j_max

    biggest_rank_that_has_a_lower_pos_than_j_in_p_terms <- -1
    for (k in 0:r_max) {
      if (k / (r_max+1) <= p_corresponding_to_j) {
        biggest_rank_that_has_a_lower_pos_than_j_in_p_terms <- k
      }else{
        break
      }
    }

    n_times <- 0

    for (k in 1:length(sortedRanks)) {
      if (sortedRanks[[k]] == biggest_rank_that_has_a_lower_pos_than_j_in_p_terms) {
        n_times = n_times + 1
      }
    }

    #since the integral needs to be 1, we need to divede the probability n_times/length(sortedRanks) by the interval length (1 / (r_max+1))
    return(n_times / length(sortedRanks) / (1 / (r_max+1)) )
}


