#' mgcvTrans: extra basis definitions for the \code{mgcv} package
#' 
#' The \code{mgcv} package offers a flexible system that allows for new basis
#' functions to be defined by users. This package introduces a handful of
#' extra basis functions, some of which are relevent for fitting smooths on
#' transformations of index variables. They can be used by specifying the
#' appropriate \code{bs} argument in one of mgcv's basis constructors (e.g.,
#' \code{s()} or \code{te()}).
#' 
#' @docType package
#' @name mgcvTrans
NULL
#> NULL