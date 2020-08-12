#' The 'CNVRG' package.
#'
#' @description This package implements Dirichlet-multinomial modeling of relative abundance data using functionality provided by the Stan software. The purpose of this package is to provide a user-friendly way to interface with Stan that is suitable for those new to R and modeling. For more regarding the modeling mathematics and computational techniques see our publication in Molecular Ecology Resources titled "Dirichlet‐multinomial modelling outperforms alternatives for analysis of microbiome and other ecological count data". DOI: 10.1111/1755-0998.13128.
#'
#' @docType package
#' @name CNVRG-package
#' @aliases CNVRG
#' @useDynLib CNVRG, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package version 2.18.2. https://mc-stan.org
#'
NULL
