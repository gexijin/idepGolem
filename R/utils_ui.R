#' ui: A fie with all functions to deal with database
#'
#'
#' @section ui functions:
#'
#'
#' 
#' @name ui
NULL


#' make_pull_down_menu
#'
#'
#' @description
#'
#'
#' @param datapath
#'
#'
#' @return
#'
#'
#' @examples
#' 
make_pull_down_menu <- function(funs) {
    return(setNames(1:length(funs), names(funs))) # nolint
}