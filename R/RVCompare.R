#' RVCompare: Compare Real Valued Random Variables
#'
#' A framework with tools to compare two random variables, and determine which of them produces lower values. It can compute the Cp and Cd of theoretical of probability distributions, as explained in E. Arza (2021) <https://github.com/EtorArza/RVCompare-paper/releases>. Given the observed samples of two random variables X_A and X_B, it can compute the confidence bands of the cumulative distributions of X'_A and X'_B (see E. Arza (2021) <https://github.com/EtorArza/RVCompare-paper> for details) based on the observed samples of X_A and X_B. Uses bootstrap and DKW-bounds to compute the confidence bands of the cumulative distributions. These two methods are described in B. Efron. (1979) <doi:10.1214/aos/1176344552> and P. Massart (1990) <doi:10.1214/aop/1176990746>.
#'
#' @docType package
#' @author Etor Arza <etorarza@gmail.com>
#' @import Rcpp ggplot2 pracma
#' @useDynLib RVCompare
#' @name RVCompare
NULL


################## Main functions ##################


#' Generate the cumulative difference-plot
#'
#' Generate the cumulative difference-plot given the observed samples using the bootstrap method.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B, it needs to have the same length as X_A.
#' @param isMinimizationProblem a boolean value where TRUE represents that lower values are preferred to larger values.
#' @param labelA (optional, default value "X_A") the label corresponding to X_A.
#' @param labelB (optional, default value "X_B") the label corresponding to X_B.
#' @param alpha (optional, default value 0.05) the error of the confidence interval. If alpha = 0.05 then we have 95 percent confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param nOfBootstrapSamples (optional, default value 1e3) how many bootstrap samples to average. Increases computation time.
#' @param ignoreMinimumLengthCheck (optional, default value FALSE) whether to skip the check for a minimum length of 100 in X_A_observed and X_A_observed.
#' @return retunrs and shows the cumulative difference-plot
#' @export
#' @examples
#'
#' ### Example 1 ###
#' X_A_observed <- rnorm(100, mean = 2, sd = 1)
#' X_B_observed <- rnorm(100, mean = 2.1, sd = 0.5)
#' \donttest{
#'  cumulative_difference_plot(X_A_observed, X_B_observed, TRUE, labelA="X_A", labelB="X_B")
#' }
#' \dontshow{
#'  cumulative_difference_plot(X_A_observed, X_B_observed, TRUE, labelA="X_A", labelB="X_B", nOfBootstrapSamples=100)
#' }
#'
#' ### Example 2 ###
#' # Comparing the optimization algorithms PL-EDA and PL-GS
#' # with 400 samples each.
#' PL_EDA_fitness <- c(
#' 52235, 52485, 52542, 52556, 52558, 52520, 52508, 52491, 52474, 52524,
#' 52414, 52428, 52413, 52457, 52437, 52449, 52534, 52531, 52476, 52434,
#' 52492, 52554, 52520, 52500, 52342, 52520, 52392, 52478, 52422, 52469,
#' 52421, 52386, 52373, 52230, 52504, 52445, 52378, 52554, 52475, 52528,
#' 52508, 52222, 52416, 52492, 52538, 52192, 52416, 52213, 52478, 52496,
#' 52444, 52524, 52501, 52495, 52415, 52151, 52440, 52390, 52428, 52438,
#' 52475, 52177, 52512, 52530, 52493, 52424, 52201, 52484, 52389, 52334,
#' 52548, 52560, 52536, 52467, 52392, 51327, 52506, 52473, 52087, 52502,
#' 52533, 52523, 52485, 52535, 52502, 52577, 52508, 52463, 52530, 52507,
#' 52472, 52400, 52511, 52528, 52532, 52526, 52421, 52442, 52532, 52505,
#' 52531, 52644, 52513, 52507, 52444, 52471, 52474, 52426, 52526, 52564,
#' 52512, 52521, 52533, 52511, 52416, 52414, 52425, 52457, 52522, 52508,
#' 52481, 52439, 52402, 52442, 52512, 52377, 52412, 52432, 52506, 52524,
#' 52488, 52494, 52531, 52471, 52616, 52482, 52499, 52386, 52492, 52484,
#' 52537, 52517, 52536, 52449, 52439, 52410, 52417, 52402, 52406, 52217,
#' 52484, 52418, 52550, 52513, 52530, 51667, 52185, 52089, 51853, 52511,
#' 52051, 52584, 52475, 52447, 52390, 52506, 52514, 52452, 52526, 52502,
#' 52422, 52411, 52171, 52437, 52323, 52488, 52546, 52505, 52563, 52457,
#' 52502, 52503, 52126, 52537, 52435, 52419, 52300, 52481, 52419, 52540,
#' 52566, 52547, 52476, 52448, 52474, 52438, 52430, 52363, 52484, 52455,
#' 52420, 52385, 52152, 52505, 52457, 52473, 52503, 52507, 52429, 52513,
#' 52433, 52538, 52416, 52479, 52501, 52485, 52429, 52395, 52503, 52195,
#' 52380, 52487, 52498, 52421, 52137, 52493, 52403, 52511, 52409, 52479,
#' 52400, 52498, 52482, 52440, 52541, 52499, 52476, 52485, 52294, 52408,
#' 52426, 52464, 52535, 52512, 52516, 52531, 52449, 52507, 52485, 52491,
#' 52499, 52414, 52403, 52398, 52548, 52536, 52410, 52549, 52454, 52534,
#' 52468, 52483, 52239, 52502, 52525, 52328, 52467, 52217, 52543, 52391,
#' 52524, 52474, 52509, 52496, 52432, 52532, 52493, 52503, 52508, 52422,
#' 52459, 52477, 52521, 52515, 52469, 52416, 52249, 52537, 52494, 52393,
#' 52057, 52513, 52452, 52458, 52518, 52520, 52524, 52531, 52439, 52530,
#' 52422, 52649, 52481, 52256, 52428, 52425, 52458, 52488, 52502, 52373,
#' 52426, 52441, 52471, 52468, 52465, 52265, 52455, 52501, 52340, 52457,
#' 52275, 52527, 52574, 52474, 52487, 52416, 52634, 52514, 52184, 52430,
#' 52462, 52392, 52529, 52178, 52495, 52438, 52539, 52430, 52459, 52312,
#' 52437, 52637, 52511, 52563, 52270, 52341, 52436, 52515, 52480, 52569,
#' 52490, 52453, 52422, 52443, 52419, 52512, 52447, 52425, 52509, 52180,
#' 52521, 52566, 52060, 52425, 52480, 52454, 52501, 52536, 52143, 52432,
#' 52451, 52548, 52508, 52561, 52515, 52502, 52468, 52373, 52511, 52516,
#' 52195, 52499, 52534, 52453, 52449, 52431, 52473, 52553, 52444, 52459,
#' 52536, 52413, 52537, 52537, 52501, 52425, 52507, 52525, 52452, 52499
#' )
#' PL_GS_fitness <- c(
#' 52476, 52211, 52493, 52484, 52499, 52500, 52476, 52483, 52431, 52483,
#' 52515, 52493, 52490, 52464, 52478, 52440, 52482, 52498, 52460, 52219,
#' 52444, 52479, 52498, 52481, 52490, 52470, 52498, 52521, 52452, 52494,
#' 52451, 52429, 52248, 52525, 52513, 52489, 52448, 52157, 52449, 52447,
#' 52476, 52535, 52464, 52453, 52493, 52438, 52489, 52462, 52219, 52223,
#' 52514, 52476, 52495, 52496, 52502, 52538, 52491, 52457, 52471, 52531,
#' 52488, 52441, 52467, 52483, 52476, 52494, 52485, 52507, 52224, 52464,
#' 52503, 52495, 52518, 52490, 52508, 52505, 52214, 52506, 52507, 52207,
#' 52531, 52492, 52515, 52497, 52476, 52490, 52436, 52495, 52437, 52494,
#' 52513, 52483, 52522, 52496, 52196, 52525, 52490, 52506, 52498, 52250,
#' 52524, 52469, 52497, 52519, 52437, 52481, 52237, 52436, 52508, 52518,
#' 52490, 52501, 52508, 52476, 52520, 52435, 52463, 52481, 52486, 52489,
#' 52482, 52496, 52499, 52443, 52497, 52464, 52514, 52476, 52498, 52496,
#' 52498, 52530, 52203, 52482, 52441, 52493, 52532, 52518, 52474, 52498,
#' 52512, 52226, 52538, 52477, 52508, 52243, 52533, 52463, 52440, 52246,
#' 52209, 52488, 52530, 52195, 52487, 52494, 52508, 52505, 52444, 52515,
#' 52499, 52428, 52498, 52244, 52520, 52463, 52187, 52484, 52517, 52504,
#' 52511, 52530, 52519, 52514, 52532, 52203, 52485, 52439, 52496, 52443,
#' 52503, 52520, 52516, 52478, 52473, 52505, 52480, 52196, 52492, 52527,
#' 52490, 52493, 52252, 52470, 52493, 52533, 52506, 52496, 52519, 52492,
#' 52509, 52530, 52213, 52499, 52492, 52528, 52499, 52526, 52521, 52488,
#' 52485, 52502, 52515, 52470, 52207, 52494, 52527, 52442, 52200, 52485,
#' 52489, 52499, 52488, 52486, 52232, 52477, 52485, 52490, 52524, 52470,
#' 52504, 52501, 52497, 52489, 52152, 52527, 52487, 52501, 52504, 52494,
#' 52484, 52213, 52449, 52490, 52525, 52476, 52540, 52463, 52200, 52471,
#' 52479, 52504, 52526, 52533, 52473, 52475, 52518, 52507, 52500, 52499,
#' 52512, 52478, 52523, 52453, 52488, 52523, 52240, 52505, 52532, 52504,
#' 52444, 52194, 52514, 52474, 52473, 52526, 52437, 52536, 52491, 52523,
#' 52529, 52535, 52453, 52522, 52519, 52446, 52500, 52490, 52459, 52467,
#' 52456, 52490, 52521, 52484, 52508, 52451, 52231, 52488, 52485, 52215,
#' 52493, 52475, 52474, 52508, 52524, 52477, 52514, 52452, 52491, 52473,
#' 52441, 52520, 52471, 52466, 52475, 52439, 52483, 52491, 52204, 52500,
#' 52488, 52489, 52519, 52495, 52448, 52453, 52466, 52462, 52489, 52471,
#' 52484, 52483, 52501, 52486, 52494, 52473, 52481, 52502, 52516, 52223,
#' 52490, 52447, 52222, 52469, 52509, 52194, 52490, 52484, 52446, 52487,
#' 52476, 52509, 52496, 52459, 52474, 52501, 52516, 52223, 52487, 52468,
#' 52534, 52522, 52474, 52227, 52450, 52506, 52193, 52429, 52496, 52493,
#' 52493, 52488, 52190, 52509, 52434, 52469, 52510, 52481, 52520, 52504,
#' 52230, 52500, 52487, 52517, 52473, 52488, 52450, 52203, 52215, 52490,
#' 52479, 52515, 52210, 52485, 52516, 52504, 52521, 52499, 52503, 52526)
#' # Considering that the LOP is a maximization problem, we need isMinimizationProblem=FALSE.
#' \donttest{
#'  cumulative_difference_plot(PL_EDA_fitness,
#'                             PL_GS_fitness,
#'                             isMinimizationProblem=FALSE,
#'                             labelA="PL-EDA",
#'                             labelB="PL-GS")
#' }
cumulative_difference_plot <- function(X_A_observed, X_B_observed, isMinimizationProblem, labelA="X_A", labelB="X_B",  alpha=0.05,  EPSILON=1e-20, nOfBootstrapSamples=1e3, ignoreMinimumLengthCheck=FALSE) {


  X_A_observed_copy <- X_A_observed
  X_B_observed_copy <- X_B_observed


  if(!isTRUE(isMinimizationProblem) && !isFALSE(isMinimizationProblem))
  {
    stop("ERROR: To compute the cumulative difference-plot, isMinimizationProblem needs to be TRUE or FALSE.
         isMinimizationProblem=TRUE needs to be used when low values are prefered to high values.
         isMinimizationProblem=FALSE needs to be used when high values are preferred to low values.")
  }
  else if (isFALSE(isMinimizationProblem))
  {
    X_A_observed_copy <- - X_A_observed_copy
    X_B_observed_copy <- - X_B_observed_copy
  }


  Y_AB_estimations <- get_Y_AB_bounds_bootstrap(X_A_observed_copy, X_B_observed_copy, alpha,  EPSILON, nOfBootstrapSamples, ignoreMinimumLengthCheck)
  res <- plot_Y_AB(Y_AB_estimations, labels = c(labelA,labelB), plotDifference = TRUE)
  return(res)
}











#' Estimate Y_A and Y_B bounds with bootstrap
#'
#' Estimate the confidence intervals for the cumulative distributions of Y_A and Y_B using bootstrap.
#' Much slower than the Dvoretzky–Kiefer–Wolfowitz approach.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B, it needs to have the same length as X_A.
#' @param alpha (optional, default value 0.05) the error of the confidence interval. If alpha = 0.05 then we have 95 percent confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param nOfBootstrapSamples (optional, default value 1e3) how many bootstrap samples to average. Increases computation time.
#' @param ignoreMinimumLengthCheck (optional, default value FALSE) wether to check for a minimum length in X_A and X_B.
#' @return Returns a list with the following fields:
#'
#' - p: values in the interval [0,1] that represent the nOfEstimationPoints points in which the densities are estimated. Useful for plotting.
#'
#' - Y_A_cumulative_estimation: an array with the estimated cumulative diustribution function of Y_A from 0 to p[[i]].
#'
#' - Y_A_cumulative_upper: an array with the upper bounds of confidence 1 - alpha of the cumulative density of Y_A
#'
#' - Y_A_cumulative_lower: an array with the lower bounds of confidence 1 - alpha of the cumulative density of Y_A
#'
#' - Y_B_cumulative_estimation: The same as Y_A_cumulative_estimation for Y_B.
#'
#' - Y_B_cumulative_upper: The same as Y_A_cumulative_upper for Y_B
#'
#' - Y_B_cumulative_lower: The same as Y_A_cumulative_lower for Y_B
#'
#' - diff_estimation: Y_A_cumulative_estimation - Y_B_cumulative_estimation
#'
#' - diff_upper: an array with the upper bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' - diff_lower: an array with the lower bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#'
#' @export
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @examples
#' library(ggplot2)
#'
#' ### Example 1 ###
#' X_A_observed <- rnorm(100, mean = 2, sd = 1)
#' X_B_observed <- rnorm(100, mean = 2.1, sd = 0.5)
#' \donttest{
#'  res <- get_Y_AB_bounds_bootstrap(X_A_observed, X_B_observed)
#' }
#' \dontshow{
#' # easier on computation for testing.
#' res <- get_Y_AB_bounds_bootstrap(X_A_observed, X_B_observed, nOfBootstrapSamples=1e2)
#' }
#' fig1 = plot_Y_AB(res, plotDifference=FALSE)+ ggplot2::ggtitle("Example 1")
#' print(fig1)
#'
#'
#' \donttest{
#' ### Example 2 ###
#' # Comparing the estimations with the actual distributions for two normal distributions.
#' ###################################
#' ## sample size = 100 ##############
#' ###################################
#' X_A_observed <- rnorm(100,mean = 1, sd = 1)
#' X_B_observed <- rnorm(100,mean = 1.3, sd = 0.5)
#' res <- get_Y_AB_bounds_bootstrap(X_A_observed, X_B_observed)
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
#' actualDistributions$Y_A_cumulative_estimation <- lm(Y_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$Y_B_cumulative_estimation <- lm(Y_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_Y_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_A_cumulative_estimation, colour = "Actual Y_A", linetype="Actual Y_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_B_cumulative_estimation, colour = "Actual Y_B", linetype="Actual Y_B")) +
#'
#' scale_colour_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="#00BFC4", "X_B"="#F8766D", "Actual Y_A"="#0000FF", "Actual Y_B"="#FF0000"))+
#'
#' scale_linetype_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="solid", "X_B"="dashed", "Actual Y_A"="solid", "Actual Y_B"="solid"))+
#'
#' ggtitle("100 samples used in the estimation")
#' print(fig)
#'
#' ###################################
#' ## sample size = 300 ##############
#' ###################################
#' X_A_observed <- rnorm(300,mean = 1, sd = 1)
#' X_B_observed <- rnorm(300,mean = 1.3, sd = 0.5)
#' res <- get_Y_AB_bounds_bootstrap(X_A_observed, X_B_observed)
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
#' actualDistributions$Y_A_cumulative_estimation <- lm(Y_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' actualDistributions$Y_B_cumulative_estimation <- lm(Y_B_cumulative_estimation ~
#'        p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'        data = actualDistributions)$fitted.values
#'
#' fig = plot_Y_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_A_cumulative_estimation, colour = "Actual Y_A", linetype="Actual Y_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_B_cumulative_estimation, colour = "Actual Y_B", linetype="Actual Y_B")) +
#'
#' scale_colour_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="#00BFC4", "X_B"="#F8766D", "Actual Y_A"="#0000FF", "Actual Y_B"="#FF0000"))+
#'
#' scale_linetype_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="solid", "X_B"="dashed", "Actual Y_A"="solid", "Actual Y_B"="solid"))+
#'
#' ggtitle("300 samples used in the estimation")
#' print(fig)
#' }
#'
get_Y_AB_bounds_bootstrap <- function(X_A_observed, X_B_observed, alpha=0.05,  EPSILON=1e-20, nOfBootstrapSamples=1e3, ignoreMinimumLengthCheck=FALSE) {

  nOfEstimationPoints = 2000

  if(sum(is.nan(c(X_A_observed, X_B_observed))) != 0)
  {
    print("ERROR: X_A_observed or X_B_observed contain NaN values.")
    stop()
  }

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
    stop()
  }

  if (alpha > 1 || alpha <= 0.0) {
    print("ERROR: alpha must be defined in the interval (0,1). It represents the error of the CI.")
    stop()
  }

  if(!ignoreMinimumLengthCheck)
  {
  if(!xHasEnoughValues(X_A_observed, 100)) {
    print("ERROR: X_A_observed does not have enough values. This means that the confidence intervals cannot be accurately computed.")
    print("At least 100 samples are required.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreMinimumLengthCheck = TRUE (not recomended!)")
    stop()
  }

  if(!xHasEnoughValues(X_B_observed, 100)) {
    print("ERROR: X_B_observed does not have enough values. This means that the confidence intervals cannot be accurately computed.")
    print("At least 100 samples are required.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreMinimumLengthCheck = TRUE (not recomended!)")
    stop()
  }
  }


  n <- length(X_A_observed)
  m <- length(X_B_observed)

  if (length(n) != length(m)) {
    print("ERROR: X_A_observed and X_B_observed need to be of equal length.")
    stop()
  }


  j_max <- nOfEstimationPoints - 1
  p <- 0:j_max / j_max # the size of the intervals is 1 / j_max.

  dataA <- matrix(0, nrow = nOfBootstrapSamples, ncol = nOfEstimationPoints)
  dataB <- matrix(0, nrow = nOfBootstrapSamples, ncol = nOfEstimationPoints)

  pb = utils::txtProgressBar(min = 1, max = nOfBootstrapSamples, initial = 1, style = 3)

  for (i in 1:nOfBootstrapSamples) {
    utils::setTxtProgressBar(pb,i)

    if(i==1)
    {
      bootStrapSampleA <- X_A_observed
      bootStrapSampleB <- X_B_observed
    }
    else
    {
      bootStrapSampleA <- sample(X_A_observed, size=n, replace=TRUE)
      bootStrapSampleB <- sample(X_B_observed, size=n, replace=TRUE)
    }
    ranksObj <- ranksOfObserved(bootStrapSampleA, bootStrapSampleB, EPSILON)

    # if(i==1 || i==2)
    # {
    #   cat("\nboot ",i, " --> ", bootStrapSampleA[1:6], "\n")
    # }

    r_max <- ranksObj$r_max
    r_max <- n

    dataA[i,] <- helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multA, j_max = j_max)
    dataB[i,] <- helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multB, j_max = j_max)
  }
  alpha_new <- 1 - sqrt(1-alpha)
  # matplot(dataA, type = c("l"),pch=1, xlab = "density dataA") #plot
  # matplot(dataB, type = c("l"),pch=1, xlab = "density dataB") #plot
  # matplot(dataA - dataB, type = c("l"),pch=1, xlab = "density dataA - dataB") #plot


  res <- list()



  res$p <- p

  # print("")
  # print("preTrapezoid")
  # print((p * 2)[1:6])
  # print(sort(dataA[1,])[1:6])
  # print(sort(dataA[2,])[1:6])


  dataA <- t(apply(dataA, 1, helperTrapezoidRule))
  dataB <- t(apply(dataB, 1, helperTrapezoidRule))

  # print("postTrapezoid")
  # print((p * 2)[1:6])
  # print(sort(dataA[1,])[1:6])
  # print(sort(dataA[2,])[1:6])
  #
  # print(which(is.na(dataA)))
  # print(which(is.na(dataB)))


  quantiles <- apply(dataA, 2, stats::quantile, probs = c(alpha_new/2, 0.5, 1.0 - alpha_new/2))



  res$Y_A_cumulative_estimation <- dataA[1,]
  res$Y_A_cumulative_lower<- quantiles[1,]
  res$Y_A_cumulative_upper <- quantiles[3,]


  quantiles <- apply(dataB, 2, stats::quantile, probs = c(alpha_new/2, 0.5, 1.0 - alpha_new/2))

  res$Y_B_cumulative_estimation <- dataB[1,]
  res$Y_B_cumulative_lower<- quantiles[1,]
  res$Y_B_cumulative_upper <- quantiles[3,]


  quantiles <- apply(dataA - dataB, 2, stats::quantile, probs = c(alpha_new/2, 0.5, 1.0 - alpha_new/2))



  res$diff_estimation <-  dataA[1,] - dataB[1,] # res$Y_A_cumulative_estimation - res$Y_B_cumulative_estimation
  res$diff_lower <- quantiles[1,] # coerceToDifferenceArea(res$p, res$Y_A_cumulative_lower - res$Y_B_cumulative_upper)
  res$diff_upper <- quantiles[3,] # coerceToDifferenceArea(res$p, res$Y_A_cumulative_upper - res$Y_B_cumulative_lower)


  return(res)
}






#' Estimate Y_A and Y_B bounds with Dvoretzky–Kiefer–Wolfowitz
#'
#'
#' Estimate the confidence intervals for the cumulative distributions of Y_A and Y_B with Dvoretzky–Kiefer–Wolfowitz.
#'
#' @param X_A_observed array of the observed samples (real values) of X_A.
#' @param X_B_observed array of the observed samples (real values) of X_B.
#' @param nOfEstimationPoints (optional, default 1000) the number of points in the interval [0,1] in which the density is estimated.
#' @param alpha (optional, default value 0.05) the error of the confidence interval. If alpha = 0.05 then we have 95 percent confidence interval.
#' @param EPSILON (optional, default value 1e-20) minimum difference between two values to be considered different.
#' @param ignoreMinimumLengthCheck (optional, default value FALSE) wether to check for a minimum length in X_A and X_B.
#' @return Returns a list with the following fields:
#'
#' - p: values in the interval [0,1] that represent the nOfEstimationPoints points in which the densities are estimated. Useful for plotting.
#'
#' - Y_A_cumulative_estimation: an array with the empirical cumulative diustribution function of Y_A from 0 to p[[i]].
#'
#' - Y_A_cumulative_upper: an array with the upper bounds of confidence 1 - alpha of the cumulative density of Y_A
#'
#' - Y_A_cumulative_lower: an array with the lower bounds of confidence 1 - alpha of the cumulative density of Y_A
#'
#' - Y_B_cumulative_estimation: The same as Y_A_cumulative_estimation for Y_B.
#'
#' - Y_B_cumulative_upper: The same as Y_A_cumulative_upper for Y_B
#'
#' - Y_B_cumulative_lower: The same as Y_A_cumulative_lower for Y_B
#'
#' - diff_estimation: Y_A_cumulative_estimation - Y_B_cumulative_estimation
#'
#' - diff_upper: an array with the upper bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' - diff_lower: an array with the lower bounds of confidence 1 - alpha of the difference between the cumulative distributions
#'
#' @export
#' @examples
#' library(ggplot2)
#' ### Example 1 ###
#' X_A_observed <- rnorm(100, mean = 2, sd = 1)
#' X_B_observed <- rnorm(100, mean = 2.1, sd = 0.5)
#' res <- get_Y_AB_bounds_DKW(X_A_observed, X_B_observed)
#' fig1 = plot_Y_AB(res, plotDifference=FALSE) + ggtitle("Example 1")
#' print(fig1)
#'
#' \donttest{
#' ### Example 2 ###
#' # Comparing the estimations with the actual distributions for two normal distributions.
#' ##################################
#' ## sample size = 100 #############
#' ##################################
#' X_A_observed <- rnorm(100,mean = 1, sd = 1)
#' X_B_observed <- rnorm(100,mean = 1.3, sd = 0.5)
#' res <- get_Y_AB_bounds_DKW(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(X_A_observed_large_sample,
#'  X_B_observed_large_sample, nOfEstimationPoints=1e4, EPSILON=1e-20)
#'
#'
#' actualDistributions$Y_A_cumulative_estimation <- lm(Y_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$Y_B_cumulative_estimation <- lm(Y_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_Y_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_A_cumulative_estimation, colour = "Actual Y_A", linetype="Actual Y_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_B_cumulative_estimation, colour = "Actual Y_B", linetype="Actual Y_B")) +
#'
#' scale_colour_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="#00BFC4", "X_B"="#F8766D", "Actual Y_A"="#0000FF", "Actual Y_B"="#FF0000"))+
#'
#' scale_linetype_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="solid", "X_B"="dashed", "Actual Y_A"="solid", "Actual Y_B"="solid"))+
#'
#' ggtitle("100 samples used in the estimation")
#' print(fig)
#'
#' ###################################
#' ## sample size = 300 ##############
#' ###################################
#' X_A_observed <- rnorm(300,mean = 1, sd = 1)
#' X_B_observed <- rnorm(300,mean = 1.3, sd = 0.5)
#' res <- get_Y_AB_bounds_DKW(X_A_observed, X_B_observed)
#'
#' X_A_observed_large_sample <- sort(rnorm(1e4, mean = 1, sd = 1))
#' X_B_observed_large_sample <- sort(rnorm(1e4, mean = 1.3, sd = 0.5))
#' actualDistributions <- getEmpiricalCumulativeDistributions(X_A_observed_large_sample,
#'  X_B_observed_large_sample, nOfEstimationPoints=1e4, EPSILON=1e-20)
#'
#'
#' actualDistributions$Y_A_cumulative_estimation <- lm(Y_A_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#' actualDistributions$Y_B_cumulative_estimation <- lm(Y_B_cumulative_estimation ~
#'         p + I(p^2) + I(p^3)+ I(p^4)+ I(p^5)+ I(p^6)+I(p^7)+ I(p^8),
#'         data = actualDistributions)$fitted.values
#'
#' fig = plot_Y_AB(res, plotDifference=FALSE) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_A_cumulative_estimation, colour = "Actual Y_A", linetype="Actual Y_A")) +
#'
#' geom_line(data=as.data.frame(actualDistributions),
#' aes(x=p, y=Y_B_cumulative_estimation, colour = "Actual Y_B", linetype="Actual Y_B")) +
#'
#' scale_colour_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="#00BFC4", "X_B"="#F8766D", "Actual Y_A"="#0000FF", "Actual Y_B"="#FF0000"))+
#'
#' scale_linetype_manual("", breaks = c("X_A", "X_B","Actual Y_A", "Actual Y_B"),
#' values = c("X_A"="solid", "X_B"="dashed", "Actual Y_A"="solid", "Actual Y_B"="solid"))+
#' ggtitle("300 samples used in the estimation")
#' print(fig)
#'}
get_Y_AB_bounds_DKW <- function(X_A_observed, X_B_observed, nOfEstimationPoints=1000, alpha=0.05,  EPSILON=1e-20, ignoreMinimumLengthCheck=FALSE) {

  if (EPSILON > 0.1 || EPSILON <= 0.0) {
    print("ERROR: EPSILON must be in the interval (0,0.1).")
    stop()

  }

  if (alpha > 1 || alpha <= 0.0) {
    print("ERROR: alpha must be defined in the interval (0,1). It represents the error of the CI.")
    stop()
  }

  if(!ignoreMinimumLengthCheck)
  {
  if(!xHasEnoughValues(X_A_observed, 100)) {
    print("ERROR: X_A_observed does not have enough values. This means that the confidence intervals cannot be accurately computed.")
    print("At least 100 samples are required.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreMinimumLengthCheck = TRUE (not recomended!)")
    stop()
  }

  if(!xHasEnoughValues(X_B_observed, 100)) {
    print("ERROR: X_B_observed does not have enough values. This means that the confidence intervals cannot be accurately computed.")
    print("At least 100 samples are required.")
    print("If you knwon what you are doing and want to proceed ignoring this error, use parameter ignoreMinimumLengthCheck = TRUE (not recomended!)")
    stop()
  }
  }

  alpha_new <- 1 - sqrt(1-alpha)

  n <- length(X_A_observed)
  bandSizeA <- sqrt( log(2 / alpha_new) / (2*n) )
  m <- length(X_B_observed)
  bandSizeB <- sqrt( log(2 / alpha_new) / (2*m) )


  if (length(n) != length(m)) {
    print("ERROR: X_A_observed and X_B_observed need to be of equal length.")
    stop()
  }

  j_max <- nOfEstimationPoints-1
  p <- 0:j_max / j_max # the size of the intervals is 1 / j_max.




  res<-getEmpiricalCumulativeDistributions(X_A_observed, X_B_observed, nOfEstimationPoints, EPSILON, trapezoid=FALSE)


  empiricalA <- res$Y_A_cumulative_estimation
  empiricalB <- res$Y_B_cumulative_estimation

  res$p <- p


  clip_to_0_1_interval <- function(y) { sapply(y, function(x) min(1,max(0,x))) }





  res$Y_A_cumulative_lower<- clip_to_0_1_interval(empiricalA - bandSizeA)
  res$Y_A_cumulative_upper <- clip_to_0_1_interval(empiricalA + bandSizeA)


  res$Y_B_cumulative_lower<- clip_to_0_1_interval(empiricalB - bandSizeB)
  res$Y_B_cumulative_upper <- clip_to_0_1_interval(empiricalB + bandSizeB)


  res$diff_estimation <- res$Y_A_cumulative_estimation - res$Y_B_cumulative_estimation
  res$diff_upper <- coerceToDifferenceArea(res$p, res$Y_A_cumulative_upper - res$Y_B_cumulative_lower)
  res$diff_lower <- coerceToDifferenceArea(res$p, res$Y_A_cumulative_lower - res$Y_B_cumulative_upper)



  return(res)


}




#' Plot the estimated cdf of Y_A and Y_B or their difference
#'
#' retunrs a ggplot2 with the estimations of  Y_A and Y_B or the difference in cumulative distribution function.
#'
#'
#' @param estimated_Y_AB_bounds the bounds estimated with \code{\link{get_Y_AB_bounds_bootstrap}} or \code{\link{get_Y_AB_bounds_DKW}}.
#' @param labels (optional, default=c("X_A","X_B")) a string vector of length 2 with the labels of X_A and X_B, in that order.
#' @param plotDifference (optional, default=TRUE) plots the difference (Y_A - Y_B) instead of each of the random variables on their own.
#' @return the ggplot figure object.
#' @export
#' @import ggplot2
#' @examples
#' ### Example 1 ###
#'
#' X_A_observed <- rnorm(800,mean = 1, sd = 1)
#' X_B_observed <- rnorm(800,mean = 1.3, sd = 0.5)
#' res <- get_Y_AB_bounds_DKW(X_A_observed, X_B_observed)
#' densitiesPlot = plot_Y_AB(res, plotDifference=TRUE)
#' print(densitiesPlot)
plot_Y_AB <- function(estimated_Y_AB_bounds, labels=c("X_A","X_B"), plotDifference=TRUE) {
  df <- data.frame(matrix(unlist(estimated_Y_AB_bounds), nrow=length(estimated_Y_AB_bounds$p), byrow=FALSE))
  colnames(df) <- names(estimated_Y_AB_bounds)

  if (length(labels) != 2) {
    print("ERROR: The length of labels shouuld be 2")
  }


  if (!is.character(labels)) {
    print("ERROR: Labels needs to be a string vector of size 2. For example, labels=c(\"Algorithm 1\", \"Algorithm 2\")")
  }

  if (labels[[1]] == labels[[2]]) {
    print("ERROR: Labels must be different from each other.")
  }


  # This annoying hack is necessary to avoid the NOTEs 'about no visible
  # binding for global variable'. This is a known problem with ggplot2, see
  # the following link:
  # https://stackoverflow.com/questions/9439256
  p <- Y_A_cumulative_estimation <- Y_B_cumulative_estimation <- NULL
  Y_A_cumulative_lower <- Y_B_cumulative_lower <- NULL
  Y_A_cumulative_upper <- Y_B_cumulative_upper <- NULL
  x <- ymin <- ymax <- NULL

  if (plotDifference) {
    diff_estimation <- estimated_Y_AB_bounds$diff_estimation
    diff_upper <- estimated_Y_AB_bounds$diff_upper
    diff_lower <- estimated_Y_AB_bounds$diff_lower
    p <- estimated_Y_AB_bounds$p
    zeros <- rep(0, length(p))
    diff_plotdf <- data.frame(p, diff_estimation, diff_lower, diff_upper,zeros)

    resPlot = ggplot2::ggplot() +
      ggplot2::xlim(0,1) +
      ggplot2::ylim(-1,1) +
      ggplot2::geom_ribbon(data = as.data.frame(list("x"=c(0.0,0.5,1.0), "ymax"=c(0.0,1.0,0.0), "ymin"=c(0.0,-1.0,0.0))), ggplot2::aes(x=x, ymin = ymin, ymax = ymax), fill = "#33ccff", alpha = 0.2) +
      ggplot2::geom_ribbon(data = diff_plotdf, ggplot2::aes(x=p, ymin = diff_lower, ymax = diff_upper), fill = "#000000", alpha = 0.2) +
      ggplot2::geom_line(data = diff_plotdf, ggplot2::aes(x=p, y=zeros), color='#53868B', size = 0.25, linetype = "dashed") +
      ggplot2::geom_line(data = diff_plotdf, ggplot2::aes(x=p, y=diff_estimation), color='black') +
      ggplot2::annotate("text", x=0.0, y=0.9, label= labels[[1]], color="#000000", hjust = 0) +
      ggplot2::annotate("text", x=0.0, y=-0.9, label= labels[[2]], color="#000000", hjust = 0) +
      ggplot2::xlab('x') +
      ggplot2::theme_minimal() +
      ggplot2::ylab('diff(x)')
    return(resPlot)

  }else{
    labelA <- labels[[1]]
    labelB <- labels[[2]]
    resPlot = ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = df, ggplot2::aes(x=p, ymin = Y_A_cumulative_lower, ymax = Y_A_cumulative_upper, fill=labelA),  alpha = 0.15) +
      ggplot2::geom_ribbon(data = df, ggplot2::aes(x=p, ymin = Y_B_cumulative_lower, ymax = Y_B_cumulative_upper, fill=labelB),  alpha = 0.15) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_A_cumulative_estimation, colour = labelA, linetype=labelA)) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_B_cumulative_estimation, colour = labelB,  linetype =labelB)) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_A_cumulative_lower, colour = labelA, linetype=labelA), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_A_cumulative_upper, colour = labelA, linetype=labelA), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_B_cumulative_lower, colour = labelB, linetype=labelB), alpha=0.4) +
      ggplot2::geom_line(data = df, ggplot2::aes(x=p, y=Y_B_cumulative_upper, colour = labelB, linetype=labelB), alpha=0.4) +
      ggplot2::scale_colour_manual("", breaks = c(labelA, labelB),  values = c("#00BFC4", "#F8766D")) +
      ggplot2::scale_linetype_manual("", breaks = c(labelA, labelB), values = c("solid", "dashed")) +
      ggplot2::scale_fill_manual("", breaks = c(labelA, labelB),  values = c("#00BFC4", "#F8766D")) +
      ggplot2::theme_minimal() +
      ggplot2::xlab('x') +
      ggplot2::ylab('cumulative probability')

    return(resPlot)
  }
}




#' Get the empirical distribution from samples.
#'
#' Given the observed sampels of X_A (or X_B) returns the empirical cumulative distribution
#' function of Y_A (or Y_B)
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
#' plot(c$p, c$Y_A_cumulative_estimation, type="l")
#' lines(x=c$p, y=c$Y_B_cumulative_estimation, col="red")
getEmpiricalCumulativeDistributions <- function(X_A_observed, X_B_observed, nOfEstimationPoints, EPSILON, trapezoid=TRUE) {



  j_max <- nOfEstimationPoints - 1
  p <- 0:j_max / j_max
  ranksObj <- ranksOfObserved(X_A_observed, X_B_observed, EPSILON)
  # X_A_ranks <- sort(ranksObj$X_A_ranks)
  # X_B_ranks <- sort(ranksObj$X_B_ranks)
  r_max <- ranksObj$r_max

  res <- list()
  if (trapezoid)
  {
    res$Y_A_cumulative_estimation <- helperTrapezoidRule(helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multA, j_max = j_max))
    res$Y_B_cumulative_estimation <- helperTrapezoidRule(helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multB, j_max = j_max))
  }
  else
  {

    res$Y_A_cumulative_estimation <- c(0, helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multA, j_max = j_max-1, cumulative = TRUE))
    res$Y_B_cumulative_estimation <- c(0, helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multB, j_max = j_max-1, cumulative = TRUE))

    # dataA <- helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multA, ranksObj$rank_interval_lengths, r_max = r_max, j_max = j_max, cumulative = FALSE)
    # dataB <- helper_from_ranks_to_integrable_values(ranksObj$rank_interval_multB, ranksObj$rank_interval_lengths, r_max = r_max, j_max = j_max, cumulative = FALSE)
    #
    #
    # print(dataA)
    # print(dataB)
    # print(r_max)
    #
    # res$Y_A_cumulative_estimation <- c(dataA)
    # res$Y_B_cumulative_estimation <- c(dataB)
    #
    # cumIntegrableValueA <- 0
    # cumIntegrableValueB <- 0
    # last_index_updated <- 1
    # index_in_data <- 1
    #
    #
    # for (i in 0:(r_max+1)) {
    #   index_in_data <- min(max(1, floor((i) / (r_max+1) * (nOfEstimationPoints)) ), length(dataA))
    #   cumIntegrableValueA <- cumIntegrableValueA + dataA[index_in_data]
    #   cumIntegrableValueB <- cumIntegrableValueB + dataB[index_in_data]
    #   while (last_index_updated <= index_in_data) {
    #     res$Y_A_cumulative_estimation[last_index_updated] = cumIntegrableValueA
    #     res$Y_B_cumulative_estimation[last_index_updated] = cumIntegrableValueB
    #     last_index_updated =  last_index_updated + 1
    #   }
    # }
    #
    #     res$Y_A_cumulative_estimation <- res$Y_A_cumulative_estimation / max(c(res$Y_A_cumulative_estimation, res$Y_B_cumulative_estimation), na.rm = TRUE)
    #     res$Y_B_cumulative_estimation <- res$Y_B_cumulative_estimation / max(c(res$Y_A_cumulative_estimation, res$Y_B_cumulative_estimation), na.rm = TRUE)
    #     res$Y_A_cumulative_estimation[0] <- 0
    #     res$Y_B_cumulative_estimation[0] <- 0
  }

  res$p <- p
  return(res)
}










################## Comparison Functions ##################


#' The dominance rate of X_A over X_B for discrete distributions, given the probability mass functions.
#'
#' Returns a real number in the interval [0,1] that represents the dominance rate of X_A over X_B.
#'
#'
#' @param pMassA The probability mass function where pMassA[[i]] is the probability of x_i, p_A(x_i).
#' @param pMassB The probability mass function where pMassB[[i]] is the probability of x_i, p_B(x_i).
#' @return Returns the dominance rate of X_A over X_B for discrete random variables.
#' @export
#' @examples
#' CdFromProbMassFunctions(c(0.2,0.6,0.2), c(0.3,0.3,0.4))
#' # > 0.6
#' # Notice how adding additional mass with the same cumulative distribution in both
#' # random variables does not change the result.
#' CdFromProbMassFunctions(c(0.2,0.6,0.2,0.2,0.2)/1.4, c(0.3,0.3,0.4,0.2,0.2)/1.4)
#' # > 0.6
CdFromProbMassFunctions <- function(pMassA, pMassB) {

  if (length(pMassA) != length(pMassB)) {
    print("ERROR: pMassA and pMassB need to be of equal length.")
    stop()
  }

  pMassA_new <- c(0, pMassA)
  pMassB_new <- c(0, pMassB)

  F_A <- cumsum(pMassA_new)
  F_B <- cumsum(pMassB_new)

  cA <- 1
  cB <- 1

  resA <- 0
  resB <- 0
  for (i in 2:length(pMassA_new)) {
    # cat("------------------\ni: ",i,"\n")
    if(F_A[i-1] == F_B[i-1] && F_A[i] == F_B[i]) # do not update resA and resB
    {
      cA = cA - pMassA_new[i]
      cB = cB - pMassB_new[i]
      next
    }
    delta_i <- NULL
    if( # F_A[i] > F_B[i] and F_A[i-1] > F_B[i-1]
      (F_A[i-1] >= F_B[i-1] && F_A[i] > F_B[i]) ||
      (F_A[i-1] > F_B[i-1]  && F_A[i] >= F_B[i])
    )
    {
      delta_i <- 0
    }else if( # F_B[i] > F_A[i] and F_B[i-1] > F_A[i-1]
      (F_B[i-1] >= F_A[i-1] && F_B[i] > F_A[i]) ||
      (F_B[i-1] > F_A[i-1]  && F_B[i] >= F_A[i])
    )
    {
      delta_i <- 1
    }
    else if(F_B[i-1] > F_A[i-1] && F_B[i] < F_A[i] ) # if F_B is higher than F_A in x_{i-1} and lower in x_{i}
    {
      delta_i <- (F_A[i-1] - F_B[i-1]) / ((F_B[i]- F_B[i-1]) - (F_A[i]- F_A[i-1]))
    }
    else if(F_B[i-1] < F_A[i-1] && F_B[i] > F_A[i] ) # if F_A is higher than F_B in x_{i-1} and lower in x_{i}
    {
      delta_i <- 1 - ((F_A[i-1] - F_B[i-1]) / ((F_B[i]- F_B[i-1]) - (F_A[i]- F_A[i-1])))
    }
    else
    {
      print("ERROR: one of the previous if's should have been true.")
      stop()
    }
    resA <- resA + pMassA_new[i]*(1 - delta_i)
    resB <- resB + pMassB_new[i]*delta_i
    # cat("delta: ",delta_i,"\n")
    # cat("resA: ",resA,"\n")
    # cat("resB: ",resB,"\n")
  }
  # cat("cA: ",cA,"\n")
  # cat("cB: ",cB,"\n")
  res <- resA / cA - resB / cB
  return(res /2 + 0.5)
}







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
    stop()
  }


  if(!isFunctionDensity(densityX_A, xlims))
  {
    print("ERROR: argument densityX_A is not a non discrete probability density function.")
    stop()
  }

  if(!isFunctionDensity(densityX_B, xlims))
  {
    print("ERROR: argument densityX_B is not a non discrete probability density function.")
    stop()
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
    stop()
  }

  if(!isFunctionDensity(densityX_A, xlims))
  {
    print("ERROR: argument densityX_A is not a non discrete probability density function.")
    stop()
  }

  if(!isFunctionDensity(densityX_B, xlims))
  {
    print("ERROR: argument densityX_B is not a non discrete probability density function.")
    stop()
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


#' Check for enough values.
#'
#' This function checks if there are at least minRequiredValues
#' values in the introduced vector.
#'
#'
#'
#' @param X the array with the values.
#' @param minRequiredValues the minimum number values required to return TRUE.
#' @return Returns TRUE if the values are OK. FALSE, if there are not enough values.
#' @export
#' @examples xHasEnoughValues(c(1,2,2,3,1,5,8,9,67,8.5,4,8.3), 6)
xHasEnoughValues <- function(X, minRequiredValues) {
  sortedX <- sort(X)

    return( length(X) >= minRequiredValues)

}







################## Probability density functions ##################

#' The probability density function of the normal distribution
#'
#' Returns the density function of the normal distribution with mean mu and standard deviation sigma.
#' The returned function is a single parameter function that returns the probability of the normal distribution in that point.
#' It is just a convinient wrapper of dnorm from the package 'stat' with some parameter checks.
#' @param mu the mean of the normal distribution.
#' @param sigma the standard deviation of the normal distribution.
#' @importFrom methods className is
#' @return Returns a callable function with a single parameter that describes the probability of the normal distribution in that point.
#' @family probability density distributions
#' @export
#' @examples
#' dist <- normalDensity(0,1)
#' dist(0)
normalDensity <- function(mu, sigma) {

  # parameter checks
  if (!is(mu, className("numeric"))) {
    print("ERROR: mu is not numeric.")
    stop()
  }
  if (!is(sigma, className("numeric"))) {
    print("ERROR: sigma is not numeric.")
    stop()
  }
  if (sigma <= 0) {
    print("ERROR: sigma needs to be positive.")
    stop()
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
    stop()
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
    stop()
  }

  n = length(densities)



  if (is.null(weights)) {
    weights = rep(1/n, n)
  }

  if (abs(sum(weights) - 1.0) > 1e-6) {
    print("ERROR: the sum of the weights must be 1.0.")
    stop()
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
    stop()
  }
  X <- 0:nIntervals / nIntervals * (xlims[[2]] - xlims[[1]]) + xlims[[1]]
  interval_size <- X[[2]] -  X[[1]]
  weights <- sapply(utils::head(X,-1) + interval_size / 2, density)
  indexes <- sample(1:nIntervals, nSamples, prob=weights, replace = TRUE)
  return (X[indexes] + stats::runif(nSamples) * interval_size)
}




################## Internal functions ##################


#' Coerce values to the square in the diffference graph
#' @param p the x values.
#' @param v the values ot be corrected.
#' @keywords internal
#' @return vector with the v values corrected
coerceToDifferenceArea <- function(p, v) {

  if(length(p) != length(v))
  {
    cat("\nERROR: the length of v should be the same of p. length(v) =",length(v), "  and length(p) =", length(p), "\n")
  }

  v_new <- v

  for (i in 1:length(p)) {
    v_new[[i]] <- min(v_new[[i]], 2*p[[i]])
    v_new[[i]] <- min(v_new[[i]], 2 - 2*p[[i]])
    v_new[[i]] <- max(v_new[[i]], -2*p[[i]])
    v_new[[i]] <- max(v_new[[i]], -2 + 2*p[[i]])
  }
  return(v_new)
}



#' Check if xlims is a tuple that represents a valid bounded interval in the real space.
#' @param xlims the tuple to be checked.
#' @importFrom methods className is
#' @return TRUE if it is a valid tuple. Otherwise prints error mesage and returns FALSE
isXlimsValid <- function(xlims) {

  if (missing(xlims)) {
    print("ERROR: xlims not set. Set the domain of definition with the parameter xlims. Example: xlims = c(-1,1) means the domain of definition is the interval (-1,1).")
    return(FALSE)
  }

  if (!is(xlims, className("numeric"))) {
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
    stop()
  }
  if (sanityChecks) {
    if(!isFunctionDensity(densityX, xlims))
    {
      print("ERROR: argument densityX_A is not a non discrete probability density function.")
      stop()
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

#' Get the ranks from the values of observed X_A and X_B. Ranks go from 0 to r_max, where r_max is the number of values in c(X_A_observed, X_B_observed)
#' @param X_A_observed array of the samples (real values) of X_A.
#' @param X_B_observed array of the samples (real values) of X_B.
#' @param EPSILON (optional, default value 1e-20) when will two values be different.
#' @return a list with three fields: X_A_ranks, X_B_ranks and r_max (the number of values minus 1).
#' @keywords internal
ranksOfObserved <- function(X_A_observed, X_B_observed, EPSILON=1e-20) {

  n = length(X_A_observed)
  m = length(X_B_observed)

  if(n != m)
  {
    print("ERROR: the length of X_A_observed and X_B_observed should be the same.")
  }

  all_values = c(X_A_observed, X_B_observed)

  order = order(all_values)
  inv_order = array(data=0, dim = length(all_values))

  for (i in 1:length(order)) {
    inv_order[order[[i]]] = i
  }


  n_repeated <- 0
  for (i in 1:(length(order)-1)) {
    if (abs(all_values[[order[[i]]]] - all_values[[order[[i+1]]]]) < EPSILON) {
      inv_order[[order[[i+1]]]] = inv_order[[order[[i]]]]
      n_repeated = n_repeated + 1
    }else{
      inv_order[[order[[i+1]]]] = inv_order[[order[[i+1]]]] - n_repeated
    }
  }


  # Modify the ranks so that the slope of $Y_A$ plus the slope of $Y_B$ is constant.

  ranksA <- inv_order[1:n]-1
  ranksB <- inv_order[n+1:m]-1

  # We compute the lowest common multiple of the repetitions
  max_rank <- max(c(ranksA, ranksB))

  lcm_on_number_of_repeated_ranks <- lowestCommonMultiple(table(c(ranksA, ranksB)))
  tabRanksAll <- table(c(ranksA, ranksB))
  tabRanksA <- table(factor(ranksA, levels=0:max_rank))
  tabRanksB <- table(factor(ranksB, levels=0:max_rank))


  # n = m
  rank_interval_multA <- c()
  rank_interval_multB <- c()
  rank_interval_lengths <- numeric(max_rank+1)

  pos_in_ranks <- 1
  for (rank in 0:max_rank) {

    rank_string <- toString(rank)

    interval_length <- tabRanksAll[[rank_string]]

    nAddNewRanksA <- tabRanksA[[rank_string]]
    nAddNewRanksB <- tabRanksB[[rank_string]]


    rank_interval_multA <- c(rank_interval_multA, rep(nAddNewRanksA, interval_length) / interval_length)
    rank_interval_multB <- c(rank_interval_multB, rep(nAddNewRanksB, interval_length) / interval_length)
    rank_interval_lengths[[rank + 1]] <- interval_length


  }


  # ranksOfObserved(c(0.1,0.2,0.5,0.5), c(0.0,0.2,0.4,0.5))
  # # Desired ouput ->
  # $rank_interval_multA
  # [1] 0 1 0.5 0.5 0 0.67 0.67 0.67
  #
  # $rank_interval_multB
  # [1] 1 0 0.5 0.5 1 0.33 0.33 0.33
  #
  # $rank_interval_lengths
  # [1] 1 1 2 1 3
  #
  # $r_max
  # [1] 5

  # the sum of rank_interval_multA and rank_interval_multB should be 1.
  # This represents that the slope is constant.
  # To achieve this, in case of repetitions, we increase the interval length
  # and the values of rank_interval_multA and rank_interval_multB represent the proportion
  # of samples in A and in B respectively.

  return(list("rank_interval_multA"=rank_interval_multA, "rank_interval_multB"=rank_interval_multB, "rank_interval_lengths"=rank_interval_lengths, "r_max"= sum(rank_interval_lengths)-1))
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
  return(as.numeric(res))
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
#' @importFrom utils tail
#' @importFrom utils head
#' @export
#' @examples
#' ### Example 1 ###
#' helperTrapezoidRule(c(1,2,3,3,3,4,5,9,3,0,1))
#' # 0.00 0.15 0.40 0.70 1.00 1.35 1.80 2.50 3.10 3.25 3.30
helperTrapezoidRule <- function(densitiesVec) {
  res <- c(0, utils::head(densitiesVec,-1) + utils::tail(densitiesVec,-1)) / 2
  return(  cumsum(res * (1 / (length(densitiesVec)-1)))  )
}





#' Helper function for get_Y_AB_bounds_bootstrap.
#'
#' The density corresponding to the position in index j is computed, given the SORTED ranks, and r_max
#' @param sortedRanks the sorted ranks of either the observed X_A or X_B.
#' @param r_max The largest rank.
#' @param j_max the largest index that will be used.
#' @param cumulative wether the integrable values should be cumulative or not. Cumulative values used in the estimation of the empirical distribution.
#' @keywords internal
#' @import Rcpp
#' @examples
#' ### Example 1 ###
#' rank_interval_mult <- c(1,0.5,1,0,0,0.5)
#' j_max <- 1000 -1
#' densities <- helper_from_ranks_to_integrable_values(rank_interval_mult, j_max, cumulative=FALSE)
#' plot(x = 0:j_max / j_max, y = densities, type="l")
#'
#' cumulative_densities <- helper_from_ranks_to_integrable_values(
#'     rank_interval_mult, j_max, cumulative=TRUE)
#' plot(x = 0:j_max / j_max, y = cumulative_densities, type="l")
#' @export
#' @return the probability density in this point
helper_from_ranks_to_integrable_values <- function(rank_interval_mult, j_max, cumulative = FALSE) {

  # since the integral needs to be 1, we need that sum_{j=0:(j_max-1)}(density in p_j * interval_length) = 1, where p_j = p[[j+1]].
  j_vec <- cpp_helper_from_ranks_to_integrable_values(rank_interval_mult, j_max)
  j_vec <- j_vec / sum(j_vec) * j_max

  if (!cumulative) {
    return(j_vec)
  }else
  {
    res <- cumsum(j_vec)
    res <- res / res[[length(res)]]
    return(res)
  }


}






