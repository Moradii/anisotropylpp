#' Nearest neighbor orientation
#'
#' Calculates the angle to the nearest neighbor point over a linear network
#'
#'
#'
#' @param X an object of class lpp
#' @param k The algorithm will find the angle towards the kth nearest neighbor
#' @param degree logical. if TRUE, output is based on angles
#' @param directed logical. if TRUE, the algorithm considers direction
#' @details this function calculates the orientation to the kth nearest neighbor for each point
#' @references Moradi, M., Mateu, J,. and Comas, C. (2020) Directional analysis for point patterns on linear networks. Stat.
#' @seealso \link[spatsat]{nndist}
#' @return a vector of angles
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' 
#' @examples
#' # generate random relaisations
#' X <- rpoislpp(0.001,branchnet)
#' nnange <- nnangle.lpp(X)


#' @import spatstat
#' @import stats
#' @export
nnangle.lpp <- function(X,k=1,degree=TRUE,directed=TRUE){
  
  ID <- nnwhich(X,k=k)
  
  ps <- psp(X$data$x,X$data$y,X$data$x[ID],X$data$y[ID],window = X$domain$window)
  ang <- angles.psp(ps,directed=directed)
  
  if(degree){
    ang[ang>0] <- (ang[ang>0]*360/(2*pi))
    ang[ang<0] <- (ang[ang<0]*360/(2*pi))+360
  }
  
  class(ang) <- "numeric"
  attr(ang,"Id") <- ID
  return(ang)
}
