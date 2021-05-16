# library(stats)
# library(pracma)


################## Main functions ##################


#' Estimate X'_A and X'_B bounds with bootstrap
#'
#' Estimate the confidence intervals for the cumulative distributions of X'_A and X'_B using bootstrap.
#' Much slower than the Dvoretzky–Kiefer–Wolfowitz approach.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B, it needs to have the same length as X_A.
#' @param nOfEstimationPoints (optional, default 100) the number of points in the interval [0,1] in which the cumulative density is estimated. Increases computation time.
#' @param alpha (optional, default value 0.2) the error of the confidence interval. If alpha = 0.05 then we have 95 percent confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param nOfBootstrapSamples (optional, default value 1e3) how many bootstrap samples to average. Increases computation time.
#' @param ignoreUniqueValuesCheck (optional, default value FALSE)
#' @return Returns a list with the following fields:
#'
#' - p: values in the interval [0,1] that represent the nOfEstimationPoints points in which the densities are estimated. Useful for plotting.
#'
#' - X_prima_A_cumulative_estimation: an array with the estimated cumulative diustribution function of X_prima_A from 0 to p[[i]].
#'
#' - X_prima_A_cumulative_upper: an array with the upper bounds of confidence 1 - alpha of the cumulative density of X_prima_A
#'
#' - X_prima_A_cumulative_lower: an array with the lower bounds of confidence 1 - alpha of the cumulative density of X_prima_A
#'
#' - X_prima_B_cumulative_estimation: The same as X_prima_A_cumulative_estimation for X'_B.
#'
#' - X_prima_B_cumulative_upper: The same as X_prima_A_cumulative_upper for X'_B
#'
#' - X_prima_B_cumulative_lower: The same as X_prima_A_cumulative_lower for X'_B
#'
#' - diff_estimation: X_prima_A_cumulative_estimation - X_prima_B_cumulative_estimation
#'
#' - diff_upper: an array with the upper bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' - diff_lower: an array with the lower bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#'
#' @export
#' @examples
#' library(ggplot2)
#'
#' ### Example 1 ###
#' X_A_observed <- c(0.13,0.21,0.13,0.11,2.2,0.12,0.5,0.14,0.21,0.17,
#'     0.11,2.0,0.12,0.50,0.14,0.16,0.2,0.23,0.6,0.11,0.18,0.113,0.1234,
#'     0.316,0.1523,0.1297,0.1123,0.139572,0.1937523)
#' X_B_observed <- c(0.71,0.12,0.19,0.17,1.5,1.0,0.5,0.41,0.11,0.16,0.01,
#'     0.31,0.34,0.64,0.14,0.13,0.09,0.21,0.29,0.36,0.41,0.13,0.142335,
#'     0.12363,0.132451,0.59217,0.157129,0.13528, 0.145)
#' \donttest{
#'  res <- get_X_prima_AB_bounds_bootstrap(X_A_observed, X_B_observed)
#' }
#' \dontshow{
#' # easier on computation for testing.
#' res <- get_X_prima_AB_bounds_bootstrap(X_A_observed, X_B_observed, nOfBootstrapSamples=1e2, nOfEstimationPoints=20)
#' }
#' fig1 = plot_X_prima_AB(res, plotDifference=FALSE)+ ggplot2::ggtitle("Example 1")
#' print(fig1)
#'
#'
#' \donttest{
#' ### Example 2 ###
#' # Comparing the estimations with the actual distributions for two normal distributions.
#' ###################################
#' ## sample size = 30 ##############
#' ###################################
#' X_A_observed <- rnorm(30,mean = 1, sd = 1)
#' X_B_observed <- rnorm(30,mean = 1.3, sd = 0.5)
#' res <- get_X_prima_AB_bounds_bootstrap(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(
#'         X_A_observed_large_sample,
#'         X_B_observed_large_sample,
#'         nOfEstimationPoints=1e4,
#'         EPSILON=1e-20)
#'
#'
#' actualDistributions$X_prima_A_cumulative_estimation <- lm(X_prima_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$X_prima_B_cumulative_estimation <- lm(X_prima_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_X_prima_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_A_cumulative_estimation, colour = "Actual X'_A", linetype="Actual X'_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_B_cumulative_estimation, colour = "Actual X'_B", linetype="Actual X'_B")) +
#'
#' scale_colour_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#' values = c("X'_A"="#F8766D", "X'_B"="#00BFC4", "Actual X'_A"="#FF0000", "Actual X'_B"="#0000FF"))+
#'
#' scale_linetype_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#' values = c("X'_A"="dashed", "X'_B"="solid", "Actual X'_A"="solid", "Actual X'_B"="solid"))+
#'
#' ggtitle("30 samples used in the estimation")
#' print(fig)
#'
#' ###################################
#' ## sample size = 300 ##############
#' ###################################
#' X_A_observed <- rnorm(300,mean = 1, sd = 1)
#' X_B_observed <- rnorm(300,mean = 1.3, sd = 0.5)
#' res <- get_X_prima_AB_bounds_bootstrap(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(
#'         X_A_observed_large_sample,
#'         X_B_observed_large_sample,
#'         nOfEstimationPoints=1e4,
#'         EPSILON=1e-20)
#'
#'
#' actualDistributions$X_prima_A_cumulative_estimation <- lm(X_prima_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' actualDistributions$X_prima_B_cumulative_estimation <- lm(X_prima_B_cumulative_estimation ~
#'        p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'        data = actualDistributions)$fitted.values
#'
#' fig = plot_X_prima_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_A_cumulative_estimation, colour = "Actual X'_A", linetype="Actual X'_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_B_cumulative_estimation, colour = "Actual X'_B", linetype="Actual X'_B")) +
#'
#' scale_colour_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#' values = c("X'_A"="#F8766D", "X'_B"="#00BFC4", "Actual X'_A"="#FF0000", "Actual X'_B"="#0000FF"))+
#'
#' scale_linetype_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#' values = c("X'_A"="dashed", "X'_B"="solid", "Actual X'_A"="solid", "Actual X'_B"="solid")) +
#'
#' ggtitle("300 samples used in the estimation")
#' print(fig)
#' }
#'
get_X_prima_AB_bounds_bootstrap <- function(X_A_observed, X_B_observed, nOfEstimationPoints=100, alpha=0.2,  EPSILON=1e-20, nOfBootstrapSamples=1e3, ignoreUniqueValuesCheck=FALSE) {

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
    return(NULL)
  }

  if (alpha > 1 || alpha <= 0.0) {
    print("ERROR: alpha must be defined in the interval (0,1). It represents the error of the CI.")
    return(NULL)
  }

  if(!xHasEnoughDiffValues(X_A_observed, EPSILON, 20)) {
    print("ERROR: X_A_observed does not have enough unique values. This means that the confidence intervals cannot be accurately computed.")
    print("Try reducing EPSILON or obtaining additional samples.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreUniqueValuesCheck = TRUE (not recomended!)")
    return(NULL)
  }

  if(!xHasEnoughDiffValues(X_B_observed, EPSILON, 20)) {
    print("ERROR: X_B_observed does not have enough unique values. This means that the confidence intervals cannot be accurately computed.")
    print("Try reducing EPSILON or obtaining additional samples.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreUniqueValuesCheck = TRUE (not recomended!)")
    return(NULL)
  }



  n <- length(X_A_observed)
  m <- length(X_B_observed)

  if (length(n) != length(m)) {
    print("ERROR: X_A_observed and X_B_observed need to be of equal length.")
    return(NULL)
  }


  nDatapointsWhereDensityEstimated <- nOfEstimationPoints - 1
  j_max <- nDatapointsWhereDensityEstimated-1
  p <- 0:j_max / j_max # the size of the intervals is 1 / j_max.

  dataA <- matrix(0, nrow = nOfBootstrapSamples, ncol = nDatapointsWhereDensityEstimated)
  dataB <- matrix(0, nrow = nOfBootstrapSamples, ncol = nDatapointsWhereDensityEstimated)

  pb = utils::txtProgressBar(min = 1, max = nOfBootstrapSamples, initial = 1, style = 3)

  for (i in 1:nOfBootstrapSamples) {
    utils::setTxtProgressBar(pb,i)
    bootStrapSampleA <- sample(X_A_observed, size= min(n,m), replace=TRUE)
    bootStrapSampleB <- sample(X_B_observed, size=min(n,m), replace=TRUE)

    ranksObj <- ranksOfObserved(bootStrapSampleA, bootStrapSampleB, EPSILON)

    X_A_ranks <- sort(ranksObj$X_A_ranks)
    X_B_ranks <- sort(ranksObj$X_B_ranks)
    r_max <- ranksObj$r_max

    dataA[i,] <- helper_from_ranks_to_integrable_values(sortedRanks = X_A_ranks, r_max = r_max, j_max = j_max)
    dataB[i,] <- helper_from_ranks_to_integrable_values(sortedRanks = X_B_ranks, r_max = r_max, j_max = j_max)
  }
  alpha_new <- 1 - sqrt(1-alpha)
  # matplot(dataA, type = c("l"),pch=1, xlab = "density dataA") #plot
  # matplot(dataB, type = c("l"),pch=1, xlab = "density dataB") #plot
  # matplot(dataA - dataB, type = c("l"),pch=1, xlab = "density dataA - dataB") #plot


  res <- list()



  res$p <- p

  dataA <- t(apply(dataA, 1, helperTrapezoidRule))
  dataB <- t(apply(dataB, 1, helperTrapezoidRule))

  quantiles <- apply(dataA, 2, stats::quantile, probs = c(alpha_new/2, 0.5, 1.0 - alpha_new/2))


  res$X_prima_A_cumulative_estimation <- quantiles[2,]
  res$X_prima_A_cumulative_lower<- quantiles[1,]
  res$X_prima_A_cumulative_upper <- quantiles[3,]


  quantiles <- apply(dataB, 2, stats::quantile, probs = c(alpha_new/2, 0.5, 1.0 - alpha_new/2))

  res$X_prima_B_cumulative_estimation <- quantiles[2,]
  res$X_prima_B_cumulative_lower<- quantiles[1,]
  res$X_prima_B_cumulative_upper <- quantiles[3,]


  res$diff_estimation <- res$X_prima_A_cumulative_estimation - res$X_prima_B_cumulative_estimation
  res$diff_upper <- res$X_prima_A_cumulative_upper - res$X_prima_B_cumulative_lower
  res$diff_lower <- res$X_prima_A_cumulative_lower - res$X_prima_B_cumulative_upper



  return(res)



}






#' Estimate X'_A and X'_B bounds with Dvoretzky–Kiefer–Wolfowitz
#'
#'
#' Estimate the confidence intervals for the cumulative distributions of X'_A and X'_B with Dvoretzky–Kiefer–Wolfowitz.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B.
#' @param nOfEstimationPoints (optional, default 1000) the number of points in the interval [0,1] in which the density is estimated.
#' @param alpha (optional, default value 0.2) the error of the confidence interval. If alpha = 0.05 then we have 95 percent confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @return Returns a list with the following fields:
#'
#' - p: values in the interval [0,1] that represent the nOfEstimationPoints points in which the densities are estimated. Useful for plotting.
#'
#' - X_prima_A_cumulative_estimation: an array with the empirical cumulative diustribution function of X_prima_A from 0 to p[[i]].
#'
#' - X_prima_A_cumulative_upper: an array with the upper bounds of confidence 1 - alpha of the cumulative density of X_prima_A
#'
#' - X_prima_A_cumulative_lower: an array with the lower bounds of confidence 1 - alpha of the cumulative density of X_prima_A
#'
#' - X_prima_B_cumulative_estimation: The same as X_prima_A_cumulative_estimation for X'_B.
#'
#' - X_prima_B_cumulative_upper: The same as X_prima_A_cumulative_upper for X'_B
#'
#' - X_prima_B_cumulative_lower: The same as X_prima_A_cumulative_lower for X'_B
#'
#' - diff_estimation: X_prima_A_cumulative_estimation - X_prima_B_cumulative_estimation
#'
#' - diff_upper: an array with the upper bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' - diff_lower: an array with the lower bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' @export
#' @examples
#' library(ggplot2)
#' ### Example 1 ###
#' X_A_observed <- c(0.13,0.21,0.13,0.11,2.2,0.12,0.5,0.14,0.21,0.17,
#'     0.11,2.0,0.12,0.50,0.14,0.16,0.2,0.23,0.6,0.11,0.18,0.113,0.1234,
#'     0.316,0.1523,0.1297,0.1123,0.139572,0.1937523)
#' X_B_observed <- c(0.71,0.12,0.19,0.17,1.5,1.0,0.5,0.41,0.11,0.16,0.01,
#'     0.31,0.34,0.64,0.14,0.13,0.09,0.21,0.29,0.36,0.41,0.13,0.142335,
#'     0.12363,0.132451,0.59217,0.157129,0.13528, 0.145)
#' res <- get_X_prima_AB_bounds_DKW(X_A_observed, X_B_observed)
#' fig1 = plot_X_prima_AB(res, plotDifference=FALSE) + ggtitle("Example 1")
#' print(fig1)
#'
#' \donttest{
#' ### Example 2 ###
#' # Comparing the estimations with the actual distributions for two normal distributions.
#' ###################################
#' ## sample size = 30 ##############
#' ###################################
#' X_A_observed <- rnorm(30,mean = 1, sd = 1)
#' X_B_observed <- rnorm(30,mean = 1.3, sd = 0.5)
#' res <- get_X_prima_AB_bounds_DKW(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(X_A_observed_large_sample,
#'  X_B_observed_large_sample, nOfEstimationPoints=1e4, EPSILON=1e-20)
#'
#'
#' actualDistributions$X_prima_A_cumulative_estimation <- lm(X_prima_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$X_prima_B_cumulative_estimation <- lm(X_prima_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_X_prima_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_A_cumulative_estimation, colour = "Actual X'_A", linetype="Actual X'_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_B_cumulative_estimation, colour = "Actual X'_B", linetype="Actual X'_B")) +
#'
#' scale_colour_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#'   values = c("X'_A"="#F8766D", "X'_B"="#00BFC4", "Actual X'_A"="#FF0000", "Actual X'_B"="#0000FF"))+
#'
#' scale_linetype_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#'  values = c("X'_A"="dashed", "X'_B"="solid", "Actual X'_A"="solid", "Actual X'_B"="solid"))+
#'
#' ggtitle("30 samples used in the estimation")
#' print(fig)
#'
#' ###################################
#' ## sample size = 300 ##############
#' ###################################
#' X_A_observed <- rnorm(300,mean = 1, sd = 1)
#' X_B_observed <- rnorm(300,mean = 1.3, sd = 0.5)
#' res <- get_X_prima_AB_bounds_DKW(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(X_A_observed_large_sample,
#'  X_B_observed_large_sample, nOfEstimationPoints=1e4, EPSILON=1e-20)
#'
#'
#' actualDistributions$X_prima_A_cumulative_estimation <- lm(X_prima_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$X_prima_B_cumulative_estimation <- lm(X_prima_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_X_prima_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_A_cumulative_estimation, colour = "Actual X'_A", linetype="Actual X'_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=X_prima_B_cumulative_estimation, colour = "Actual X'_B", linetype="Actual X'_B")) +
#'
#' scale_colour_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#'   values = c("X'_A"="#F8766D", "X'_B"="#00BFC4", "Actual X'_A"="#FF0000", "Actual X'_B"="#0000FF"))+
#'
#'
#' scale_linetype_manual("", breaks = c("X'_A", "X'_B","Actual X'_A", "Actual X'_B"),
#'  values = c("X'_A"="dashed", "X'_B"="solid", "Actual X'_A"="solid", "Actual X'_B"="solid")) +
#' ggtitle("300 samples used in the estimation")
#' print(fig)
#'}
get_X_prima_AB_bounds_DKW <- function(X_A_observed, X_B_observed, nOfEstimationPoints=1000, alpha=0.2,  EPSILON=1e-20) {

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
    return(NULL)

  }

  if (alpha > 1 || alpha <= 0.0) {
    print("ERROR: alpha must be defined in the interval (0,1). It represents the error of the CI.")
    return(NULL)
  }


  if(!xHasEnoughDiffValues(X_A_observed, EPSILON, 20)) {
    print("ERROR: X_A_observed does not have enough unique values. This means that the confidence intervals cannot be accurately computed.")
    print("Try reducing EPSILON or obtaining additional samples.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreUniqueValuesCheck = TRUE (not recomended!)")
    return(NULL)
  }

  if(!xHasEnoughDiffValues(X_B_observed, EPSILON, 20)) {
    print("ERROR: X_B_observed does not have enough unique values. This means that the confidence intervals cannot be accurately computed.")
    print("Try reducing EPSILON or obtaining additional samples.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreUniqueValuesCheck = TRUE (not recomended!)")
    return(NULL)
}


  alpha_new <- 1 - sqrt(1-alpha)

  n <- length(X_A_observed)
  bandSizeA <- sqrt( log(2 / alpha_new) / (2*n) )
  m <- length(X_B_observed)
  bandSizeB <- sqrt( log(2 / alpha_new) / (2*m) )


  if (length(n) != length(m)) {
    print("ERROR: X_A_observed and X_B_observed need to be of equal length.")
    return(NULL)
  }

  nDatapointsWhereDensityEstimated <- nOfEstimationPoints - 1
  j_max <- nDatapointsWhereDensityEstimated-1
  p <- 0:j_max / j_max # the size of the intervals is 1 / j_max.


  ranksObj <- ranksOfObserved(X_A_observed, X_B_observed, EPSILON)

  X_A_ranks <- sort(ranksObj$X_A_ranks)
  X_B_ranks <- sort(ranksObj$X_B_ranks)
  r_max <- ranksObj$r_max

  res<-getEmpiricalCumulativeDistributions(X_A_observed, X_B_observed, nOfEstimationPoints, EPSILON, trapezoid=FALSE)


  empiricalA <- res$X_prima_A_cumulative_estimation
  empiricalB <- res$X_prima_B_cumulative_estimation

  res$p <- p


  clip_to_0_1_interval <- function(y) { sapply(y, function(x) min(1,max(0,x))) }





  res$X_prima_A_cumulative_lower<- clip_to_0_1_interval(empiricalA - bandSizeA)
  res$X_prima_A_cumulative_upper <- clip_to_0_1_interval(empiricalA + bandSizeA)


  res$X_prima_B_cumulative_lower<- clip_to_0_1_interval(empiricalB - bandSizeB)
  res$X_prima_B_cumulative_upper <- clip_to_0_1_interval(empiricalB + bandSizeB)


  res$diff_estimation <- res$X_prima_A_cumulative_estimation - res$X_prima_B_cumulative_estimation
  res$diff_upper <- res$X_prima_A_cumulative_upper - res$X_prima_B_cumulative_lower
  res$diff_lower <- res$X_prima_A_cumulative_lower - res$X_prima_B_cumulative_upper



  return(res)


}




#' Plot the estimated cdf of X'_A and X'_B or their difference
#'
#' retunrs a ggplot2 with the estimations of  X'_A and X'_B or the difference in cumulative distribution function.
#'
#'
#' @param estimated_X_prima_AB_bounds the bounds estimated with \code{\link{get_X_prima_AB_bounds_bootstrap}} or \code{\link{get_X_prima_AB_bounds_DKW}}.
#' @param labels (optional, c("X'_A","X'_B")) a string vector of length 2 with the labels of X_A and X_B, in that order.
#' @param plotDifference (optional, default=TRUE) plots the difference (X'_A - X'_B) instead of each of the random variables on their own.
#' @return the ggplot figure object.
#' @export
#' @import ggplot2
#' @examples
#' ### Example 1 ###
#'
#' X_A_observed <- rnorm(800,mean = 1, sd = 1)
#' X_B_observed <- rnorm(800,mean = 1.3, sd = 0.5)
#' res <- get_X_prima_AB_bounds_DKW(X_A_observed, X_B_observed)
#' densitiesPlot = plot_X_prima_AB(res, plotDifference=TRUE)
#' print(densitiesPlot)
plot_X_prima_AB <- function(estimated_X_prima_AB_bounds, labels=c("X'_A","X'_B"), plotDifference=TRUE) {
  df <- data.frame(matrix(unlist(estimated_X_prima_AB_bounds), nrow=length(estimated_X_prima_AB_bounds$p), byrow=FALSE))
  colnames(df) <- names(estimated_X_prima_AB_bounds)

  if (length(labels) != 2) {
    print("ERROR: The length of labels shouuld be 2")
  }

  if (class(labels) != "character") {
    print("ERROR: Labels needs to be a string vector of size 2. For example, labels=c(\"Algorithm 1\", \"Algorithm 2\")")
  }

  if (labels[[1]] == labels[[2]]) {
    print("ERROR: Labels must be different from each other.")
  }


  # This annoying hack is necessary to avoid the NOTEs 'about no visible
  # binding for global variable'. This is a known problem with ggplot2, see
  # the following link:
  # https://stackoverflow.com/questions/9439256
  p <- X_prima_A_cumulative_estimation <- X_prima_B_cumulative_estimation <- NULL
  X_prima_A_cumulative_lower <- X_prima_B_cumulative_lower <- NULL
  X_prima_A_cumulative_upper <- X_prima_B_cumulative_upper <- NULL
  x <- ymin <- ymax <- NULL

  if (plotDifference) {
    diff_estimation <- estimated_X_prima_AB_bounds$X_prima_A_cumulative_estimation - estimated_X_prima_AB_bounds$X_prima_B_cumulative_estimation
    diff_upper <- estimated_X_prima_AB_bounds$X_prima_A_cumulative_upper - estimated_X_prima_AB_bounds$X_prima_B_cumulative_lower
    diff_lower <- estimated_X_prima_AB_bounds$X_prima_A_cumulative_lower - estimated_X_prima_AB_bounds$X_prima_B_cumulative_upper
    p <- estimated_X_prima_AB_bounds$p
    diff_plotdf <- data.frame(p, diff_estimation, diff_lower, diff_upper)

    resPlot = ggplot2::ggplot() +
      ggplot2::xlim(0,1) +
      ggplot2::ylim(-1,1) +
      ggplot2::geom_ribbon(data = as.data.frame(list("x"=c(0.0,0.5,1.0), "ymax"=c(0.0,1.0,0.0), "ymin"=c(0.0,-1.0,0.0))), ggplot2::aes(x=x, ymin = ymin, ymax = ymax), fill = "#000000", alpha = 0.075) +
      ggplot2::geom_ribbon(data = diff_plotdf, ggplot2::aes(x=p, ymin = diff_lower, ymax = diff_upper), fill = "#000000", alpha = 0.25) +
      ggplot2::geom_line(data = diff_plotdf, ggplot2::aes(x=p, y=diff_estimation), color='black') +
      ggplot2::annotate("text", x=0.0, y=0.9, label= labels[[1]], color="#000000", hjust = 0) +
      ggplot2::annotate("text", x=0.0, y=-0.9, label= labels[[2]], color="#000000", hjust = 0) +
      ggplot2::xlab('x') +
      ggplot2::theme_minimal() +
      ggplot2::ylab('difference in cumulative probability')
    return(resPlot)

  }else{
    labelA <- labels[[1]]
    labelB <- labels[[2]]
    resPlot = ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = df, ggplot2::aes(x=p, ymin = X_prima_A_cumulative_lower, ymax = X_prima_A_cumulative_upper), fill="#00BFC4",  alpha = 0.15) +
      ggplot2::geom_ribbon(data = df, ggplot2::aes(x=p, ymin = X_prima_B_cumulative_lower, ymax = X_prima_B_cumulative_upper), fill="#F8766D",  alpha = 0.15) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_A_cumulative_estimation, colour = labelA, linetype=labelA)) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_B_cumulative_estimation, colour = labelB,  linetype =labelB)) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_A_cumulative_lower, colour = labelA, linetype=labelA), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_A_cumulative_upper, colour = labelA, linetype=labelA), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_B_cumulative_lower, colour = labelB, linetype=labelB), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=X_prima_B_cumulative_upper, colour = labelB, linetype=labelB), alpha=0.4) +
      ggplot2::scale_colour_manual("", breaks = c(labelA, labelB),  values = c("#00BFC4", "#F8766D")) +
      ggplot2::scale_linetype_manual("", breaks = c(labelA, labelB), values = c("solid", "dashed")) +
      ggplot2::xlab('x') +
      ggplot2::ylab('cumulative probability')

    return(resPlot)
  }
}



#' Get the empirical distribution from samples.
#'
#' Given the observed sampels of X_A (or X_B) returns the empirical cumulative distribution
#' function of X'_A (or X'_B)
#'
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B.
#' @param nOfEstimationPoints the number of points in the interval [0,1] in which the cumulative density is estimated + 2.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param trapezoid (optional, default TRUE) if trapezoid=FALSE the non smooth empirical distribution is given. This is
#' what the WDK uses the empirical as the estimation.
#' @return a list with two fields: the empirical distributions of X'A and X'B.
#' @export
#' @examples
#' ### Example 1 ###
#' c <- getEmpiricalCumulativeDistributions(c(1:5),c(1:3,2:3), 170, EPSILON=1e-20, trapezoid=FALSE)
#' plot(c$p, c$X_prima_A_cumulative_estimation, type="l")
#' lines(x=c$p, y=c$X_prima_B_cumulative_estimation, col="red")
getEmpiricalCumulativeDistributions <- function(X_A_observed, X_B_observed, nOfEstimationPoints, EPSILON, trapezoid=TRUE) {

  j_max <- nOfEstimationPoints -2
  ranksObj <- ranksOfObserved(X_A_observed, X_B_observed, EPSILON)
  X_A_ranks <- sort(ranksObj$X_A_ranks)
  X_B_ranks <- sort(ranksObj$X_B_ranks)
  r_max <- ranksObj$r_max

  res <- list()
  if (trapezoid)
  {
    res$X_prima_A_cumulative_estimation <- helperTrapezoidRule(helper_from_ranks_to_integrable_values(sortedRanks = X_A_ranks, r_max = r_max, j_max = j_max))
    res$X_prima_B_cumulative_estimation <- helperTrapezoidRule(helper_from_ranks_to_integrable_values(sortedRanks = X_B_ranks, r_max = r_max, j_max = j_max))
  }
  else
  {
  #for A
    positions <- floor(X_A_ranks / (r_max) * j_max)
    res$X_prima_A_cumulative_estimation <- array(0, dim=j_max + 1)
    tabA <- table(X_A_ranks)



    for (rank_idx in 1:length(X_A_ranks)) {
      rank <- X_A_ranks[[rank_idx]]
      nreps <- as.numeric(tabA[ names(tabA)==rank ])

      if (length(nreps)==0) {
        nreps <- 0
      }
      res$X_prima_A_cumulative_estimation[[    positions[[ rank_idx ]]  + 1   ]] <- nreps
    }
    res$X_prima_A_cumulative_estimation <- cumsum(res$X_prima_A_cumulative_estimation) / sum(res$X_prima_A_cumulative_estimation)
    res$X_prima_A_cumulative_estimation[[1]] <- 0




    # for B
    positions <- floor(X_B_ranks / (r_max) * j_max)
    res$X_prima_B_cumulative_estimation <- array(0, dim=j_max + 1)
    tabB <- table(X_B_ranks)


    for (rank_idx in 1:length(X_B_ranks)) {
      rank <- X_B_ranks[[rank_idx]]
      nreps <- as.numeric(tabB[ names(tabB)==rank ])

      if (length(nreps)==0) {
        nreps <- 0
      }
      res$X_prima_B_cumulative_estimation[[    positions[[ rank_idx ]]  + 1   ]] <- nreps
    }
    res$X_prima_B_cumulative_estimation <- cumsum(res$X_prima_B_cumulative_estimation) / sum(res$X_prima_B_cumulative_estimation)
    res$X_prima_B_cumulative_estimation[[1]] <- 0
  }



  res$p <- 0:j_max / j_max
  return(res)
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
#' @param EPSILON (optional, default = 1e-3) minimum difference between two values.
#' @import pracma
#' @return Returns the dominance rate of X_A over X_B.
#' @seealso \code{\link{CpFromDensities}}
#' @export
#' @examples
#' \donttest{
### Example 1 ###
#' # If two symmetric distributions are centered in the same point (x = 0 in
#' # this case), then their Cd will be 0.5.
#'
#' densityX_A <- normalDensity(0,1)
#' densityX_B <- uniformDensity(c(-2,2))
#' CdFromDensities(densityX_A, densityX_B, c(-5,5))
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
#' densityX_B <- mixtureDensity(c(normalDensity(0.05025,0.0015),
#'                                normalDensity(0.04525, 0.0015)),
#'                                weights = c(1 - tau, tau))
#' plot(densityX_A, from=0.03, to=0.07, type="l",  col="red", xlab="x", ylab="probability density")
#' curve(densityX_B, add=TRUE, col="blue", type="l", lty=2)
#' Cd <- CdFromDensities(densityX_A, densityX_B, c(.03,.07))
#' mtext(paste("Cd(X_A, X_B) =", format(round(Cd, 3), nsmall = 3)), side=3) # add Cd to plot as text
#' legend(x = c(0.0325, 0.045), y = c(200, 250),legend=c("X_A", "X_B"),
#'                                              col=c("red", "blue"),
#'                                              lty=1:2,
#'                                              cex=0.8) # add legend
#'
#'
#' ### Example 4 ###
#' # The dominance factor ignores the mass of the probability where the
#' # distribution functinos are equal.
#' densityX_A <- uniformDensity(c(0.1, 0.3))
#' densityX_B <- uniformDensity(c(-0.2,0.5))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-2,2))
#'
#' densityX_A <- mixtureDensity(c(uniformDensity(c(0.1,0.3)), uniformDensity(c(-1,-0.5))))
#' densityX_B <- mixtureDensity(c(uniformDensity(c(-0.2,0.5)), uniformDensity(c(-1,-0.5))))
#' CdFromDensities(densityX_A, densityX_B, xlims = c(-2,2))
#' }
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
  cA = pracma::integral(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_A(x)}, xmin=xlims[[1]], xmax=xlims[[2]], method = "Simpson") # the cA in the paper is cA^-1
  cB = pracma::integral(function(x) {as.integer(abs(cumX_A(x) - cumX_B(x)) > EPSILON) * densityX_B(x)}, xmin=xlims[[1]], xmax=xlims[[2]], method = "Simpson") # the cB in the paper is cB^-1

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

  return(0.5 *pracma::integral(f_to_integrate, xmin = xlims[[1]], xmax = xlims[[2]], method = "Simpson") + 0.5)

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
#' # If two symmetric distributions are centered in the same point (x = 0 in
#' # this case), then their Cp will be 0.5.
#' densityX_A <- normalDensity(0,1)
#' densityX_B <- uniformDensity(c(-2,2))
#' Cp = CpFromDensities(densityX_A, densityX_B, c(-5,5))
#' plot(densityX_A, from=-5, to=5, type="l",  col="red", xlab="x", ylab="probability density")
#' curve(densityX_B, add=TRUE, col="blue", type="l", lty=2)
#' mtext(paste("Cp(X_A, X_B) =", format(round(Cp, 3), nsmall = 3)), side=3) # add Cp to plot as text
#' legend(x = c(-4.5, -2), y = c(0.325, 0.4),legend=c("X_A", "X_B"),
#'                                           col=c("red", "blue"),
#'                                           lty=1:2, cex=0.8) # add legend
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

  f_to_integrate = function(y) { sapply(y, function(x) {densityX_A(x) * stats::integrate(densityX_B, lower=x, upper=xlims[[2]])$value}) }

  return(stats::integrate(f_to_integrate, lower = xlims[[1]], upper = xlims[[2]])$value)
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
#' @import pracma
#' @export
#' @examples
#' dist1 <- normalDensity(0,1)
#' # the integral of the density of the normal distribution is too low in the interval (-2,2)
#' isFunctionDensity(dist1, c(-2,2))
#' isFunctionDensity(dist1, c(-5,5)) # it is close enough from 1 in the interval (-5,5)

#' dist2 <- uniformDensity(c(0,1))
#' isFunctionDensity(dist2, xlims=c(-2,2))
#' isFunctionDensity(dist2, xlims=c(0.5,2)) # the integral is not 1
#'
#' dist3 <- function(x) 0.5/sqrt(x)
#' # The integral of the function being 1 is not enough to be considered a density function.
#' # It also needs to be boounded.
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


  integrand = pracma::integral(f, xmin = xlims[[1]], xmax = xlims[[2]])
  if (  !(abs(integrand  -1 ) < tol)  )  {
    print(paste("ERROR: the integral of a density function in its domain must be 1. The value of the integral was ", toString(integrand), " instead.", sep=""))
    return(FALSE)
  }

  return(TRUE)
}


#' Check for enough unique values.
#'
#' This function checks if there are at least minRequiredDiffValues unique
#' values in the introduced vector.
#'
#'
#'
#' @param X the array with the values.
#' @param EPSILON when will two values be considered different.
#' @param minRequiredDiffValues the minimum number of different values required to return TRUE.
#' @return Returns TRUE if the values are OK. FALSE, if there are not enough unique values.
#' @export
#' @examples xHasEnoughDiffValues(c(1,2,2,3,1,5,8,9,67,8.5,4,8.3), 1e-9, 6)
xHasEnoughDiffValues <- function(X, EPSILON, minRequiredDiffValues) {
  sortedX <- sort(X)

  nDiff <- 0
  for (i in 1:(length(sortedX)-1)) {
    if( abs(sortedX[[i]] - sortedX[[i+1]]) > EPSILON )
    {
      nDiff = nDiff + 1
    }
  }

  return( nDiff >= minRequiredDiffValues)

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

  return(function(x) stats::dnorm(x, mean = mu, sd = sigma))
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

  return(function(x) stats::dunif(x, min=xlims[[1]], max=xlims[[2]]))
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
#' plot(dist1, xlim = c(-5,5), xlab="x", ylab = "Probability density",
#'      main="Mixture of two Gaussians with equal weights", cex.main=0.85)
#'
#' dist2 <- mixtureDensity(c(normalDensity(-2,1), normalDensity(2,1)), weights=c(0.8,0.2))
#' plot(dist2, xlim = c(-5,5), xlab="x", ylab = "Probability density",
#'      main="Mixture of two Gaussians with different weights", cex.main=0.85)
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




#' Get sample given the density function
#'
#' Returns an array with samples given the probability density function.
#' @param density the probability density function.
#' @param nSamples the number of samples to generate.
#' @param xlims the domain of definition of the random variable.
#' @param nIntervals (optional, default = 1e4) the number of intervals from which to draw samples. A higher value implies more accuracy but also more computation time.
#' @return Returns an array of samples.
#' @export
#' @examples
#' normDens <- normalDensity(0,1)
#' samples <- sampleFromDensity(normDens, 1e4, c(-4,4))
#' hist(samples,  breaks=20)
sampleFromDensity <- function(density, nSamples, xlims, nIntervals=1e5) {

  # parameter checks
  if (!isXlimsValid(xlims)) {
    print("ERROR: xlims in not correctly defined.")
    return(NULL)
  }
  X <- 0:nIntervals / nIntervals * (xlims[[2]] - xlims[[1]]) + xlims[[1]]
  interval_size <- X[[2]] -  X[[1]]
  weights <- sapply(utils::head(X,-1) + interval_size / 2, density)
  indexes <- sample(1:nIntervals, nSamples, prob=weights, replace = TRUE)
  return (X[indexes] + stats::runif(nSamples) * interval_size)
}




################## Internal functions ##################

#' Check if xlims is a tuple that represents a valid bounded interval in the real space.
#' @param xlims the tuple to be checked.
#' @return TRUE if it is a valid tuple. Otherwise prints error mesage and returns FALSE
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
#' @import pracma
#' @return a callable function representing the cumulative distribution.
#' @keywords internal
#' @export
#' @examples
#' cumulativeProbability <- cumulativeFromDensity(normalDensity(0,1), c(-4,4), FALSE)
#' x <- seq(-4, 4, length.out=101)
#' plot(x, cumulativeProbability(x), type = "l")
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
    return(pracma::integral(densityX, xmin=xlims[[1]], xmax=x, method="Simpson"))
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


  # Modify the ranks so that the slope of $X'_A$ plus the slope of $X'_B$ is constant.

  ranksA <- inv_order[1:n]-1
  ranksB <- inv_order[n+1:m]-1

  # We compute the lowest common multiple of the repetitions
  max_rank <- max(c(ranksA, ranksB))

  lcm_on_number_of_repeated_ranks <- lowestCommonMultiple(table(c(ranksA, ranksB)))
  tabRanksAll <- table(c(ranksA, ranksB))
  tabRanksA <- table(factor(ranksA, levels=0:max_rank))
  tabRanksB <- table(factor(ranksB, levels=0:max_rank))

  newRanksA <- c()
  newRanksB <- c()

  current_rank <- 0
  for (rank in 0:max_rank) {
    interval_length <- tabRanksAll[toString(rank)]

    nAddNewRanksA <- lcm_on_number_of_repeated_ranks / tabRanksAll[toString(rank)] * tabRanksA[toString(rank)]
    nAddNewRanksB <- lcm_on_number_of_repeated_ranks / tabRanksAll[toString(rank)] * tabRanksB[toString(rank)]

    newRanksA <- c(newRanksA, sort(rep(current_rank:(current_rank+interval_length-1), nAddNewRanksA)))
    newRanksB <- c(newRanksB, sort(rep(current_rank:(current_rank+interval_length-1), nAddNewRanksB)))
    current_rank = current_rank + interval_length
  }




  # ranksOfObserved(c(0.1,0.2,0.5), c(0.2,0.4,0.5))
  # # Desired ouput ->
  # $X_A_ranks
  # [1] 0 0  1 2       4 5
  # $X_B_ranks
  # [1]      1 2  3 3  4 5

  # # explanation -> we need the slope two be half when two ranks are shared in both positions.
  # thus, we duplicate the slope in the rest of the positions. What about 2 reps in A and 3 in B shared?

  return(list("X_A_ranks"=newRanksA, "X_B_ranks"=newRanksB, "r_max"= max(c(newRanksA,newRanksB))))
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



#' Get the lowest common divisor
#' @param integerArray an array of integers.
#' @return a integer that is the lcd of the numbers in integerArray
#' @keywords internal
lowestCommonMultiple <- function(integerArray) {
  res <-1
  uniueIntegerArray <- unique(integerArray)
  for (value in uniueIntegerArray) {
    res <- pracma::Lcm(res, value)
  }
  return(res)
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
#' matplot(x,cbind(res$X_prima_A(x),res$X_prima_B(x)),type="l",col=c("red","blue"),
#'                     ylab='Probability density')
#' legend(x = c(0.7, 1.0), y = c(2.0, 2.5),legend=c("X'_A", "X'_B"), col=c("red", "blue"),
#'                     lty=1:2, cex=0.8) # add legend
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




#' Helper function to compute the integrals in each interval.
#'
#' Computes a vector in which each index i is the integral in the interval (0, p[[i]])
#' of the function described by the densityVec
#' Uses the trapezoidal rule # https://en.wikipedia.org/wiki/Trapezoidal_rule
#' to integrate the values in the interval [0,1]. The x corexponding to the
#' values (the f(x)) are assumed to be equidistantly distributed in the interval,
#' where the x corresponding to densitiesVec[[1]] is located in 0.0 andthe x
#' corresponding to densitiesVec[[length(densitiesVec)]] is located in 1.0
#'
#' @param densitiesVec the vector of values to be integrated
#' @keywords internal
#' @return a vector in which each index i is the integral in the interval
#' (0, p[[i]]). Consequently, the first element in the vector returned
#' will be 0, since p[[1]] = 0 does not exist.
#' @export
#' @examples
#' ### Example 1 ###
#' helperTrapezoidRule(c(1,2,3,3,3,4,5,9,3,0,1))
#' # 0.00 0.15 0.40 0.70 1.00 1.35 1.80 2.50 3.10 3.25 3.30
helperTrapezoidRule <- function(densitiesVec) {
  res <- c(0, utils::head(densitiesVec,-1) + utils::tail(densitiesVec,-1)) / 2
  return(  cumsum(res * (1 / (length(densitiesVec)-1)))  )
}





#' Helper function for get_X_prima_AB_bounds_bootstrap.
#'
#' The density corresponding to the position in index j is computed, given the SORTED ranks, and r_max
#' @param sortedRanks the sorted ranks of either the observed X_A or X_B.
#' @param r_max The largest rank.
#' @param j_max the largest index that will be used. its value is nDatapointsWhereDensityEstimated-1
#' @keywords internal
#' @examples
#' ### Example 1 ###
#' j_max <- 12
#' r_max <- 6
#' sortedRanks <- c(0,0,0,0,1,1,1,1,1,1,1,3,4)
#' densities <- helper_from_ranks_to_integrable_values(
#'              sortedRanks=sortedRanks, r_max=r_max, j_max=j_max)
#' plot(x = 0:j_max / j_max, y = densities, type="l")
#' # 0.9347826 0.9782609 1.0000000 1.0000000 1.0000000 1.0000000
#' print(utils::tail(helperTrapezoidRule(densities)))
#' plot(x = 0:j_max / j_max, y = helperTrapezoidRule(densities), type="l")
#'
#' ### Example 2 ###
#' j_max <- 12
#' r_max <- 19
#' sortedRanks <- c(0,0,1,1,3,5,6,18,19)
#' densities <- helper_from_ranks_to_integrable_values(
#'                                  sortedRanks=sortedRanks, r_max=r_max, j_max=j_max)
#' plot(x = 0:j_max / j_max, y = densities, type="l")
#' # 0.8000000 0.8000000 0.8000000 0.8000000 0.8666667 1.0000000
#' print(utils::tail(helperTrapezoidRule(densities)))
#' plot(x = 0:j_max / j_max, y = helperTrapezoidRule(densities), type="l")
#'
#' ### Example 3 ###
#' j_max <- 12
#' r_max <- 8
#' sortedRanks <- c(1,1,3,5,6)
#' densities <- helper_from_ranks_to_integrable_values(
#'              sortedRanks=sortedRanks, r_max=r_max, j_max=j_max)
#' plot(x = 0:j_max / j_max, y = densities, type="l")
#' # 0.6428571 0.7857143 0.9285714 1.0000000 1.0000000 1.0000000
#' print(utils::tail(helperTrapezoidRule(densities)))
#' plot(x = 0:j_max / j_max, y = helperTrapezoidRule(densities), type="l")
#' @export
#' @return the probability density in this point
helper_from_ranks_to_integrable_values <- function(sortedRanks, r_max, j_max) {


  j_vec = array(0, dim=j_max+1)

  last_biggest_index <- 0
  n_times_last <- 0
  last_k_index <- 1
  # j goes from 0 to j_max
  for (j in 0:j_max) {
    p_corresponding_to_j <- j / j_max

    biggest_rank_that_has_a_lower_pos_than_j_in_p_terms <- -1
    k_was_updated <- FALSE
    for (k in last_biggest_index:r_max) {
      if (k / (r_max+1) <= p_corresponding_to_j) {
        biggest_rank_that_has_a_lower_pos_than_j_in_p_terms <- k
        k_was_updated <- TRUE
      }else{
        last_biggest_index <- biggest_rank_that_has_a_lower_pos_than_j_in_p_terms
        break
      }
    }

    n_times <- 0

    if (!k_was_updated) {
      n_times <- n_times_last
    }else{
      first_k_index_not_found <- TRUE
      for (k_index in last_k_index:length(sortedRanks)) {
        if (sortedRanks[[k_index]] == biggest_rank_that_has_a_lower_pos_than_j_in_p_terms) {
          if(first_k_index_not_found)
          {
            last_k_index <- k_index
            first_k_index_not_found <- FALSE
          }
          n_times = n_times + 1
        }else if(!first_k_index_not_found){
          break
        }
      }
    }
    j_vec[[j+1]] <- n_times
    n_times_last <- n_times
  }

  #since the integral needs to be 1, we need that sum_{j=0:(j_max-1)}(density in p_j * interval_length) = 1, where p_j = p[[j+1]].
  return(j_vec / utils::tail(helperTrapezoidRule(j_vec),1))
}




