#' RVCompare: Compare Real Valued Random Variables
#'
#' A framework with tools to compare two random variables, and determine which of them produces lower values. It can compute the Cp and Cd of theoretical of probability distributions, as explained in E. Arza (2021) <https://github.com/EtorArza/RVCompare-paper/releases>. Given the observed samples of two random variables X_A and X_B, it can compute the confidence bands of the cumulative distributions of X'_A and X'_B (see E. Arza (2021) <https://github.com/EtorArza/RVCompare-paper> for details) based on the observed samples of X_A and X_B. Uses bootstrap and DKW-bounds to compute the confidence bands of the cumulative distributions. These two methods are described in B. Efron. (1979) <doi:10.1214/aos/1176344552> and P. Massart (1990) <doi:10.1214/aop/1176990746>.
#'
#' @docType package
#' @author Etor Arza <etorarza@gmail.com>
#' @import Rcpp ggplot2 pracma
#' @importFrom Rcpp evalCpp
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom utils head
#' @importFrom utils tail
#' @useDynLib RVCompare
#' @name RVCompare
NULL
