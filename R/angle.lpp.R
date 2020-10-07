#' Orientation angles of lines connecting points over a linear network
#'
#' Calculates the orientation angle of lines connecting points over a linear network
#'
#'
#'
#' @param X an object of class lpp
#' @param degree logical. if TRUE, output is based on angles from 0 to 360
#' @param directed logical. if TRUE, the algorithm respects the sense of direction, see \link[spatstat]{angles.psp}
#' 
#' @details This function calculates the pairwise orientation between points over a linear network.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' @seealso \link[spatstat]{angles.psp}
#' @return a matrix of angles
#' @references Moradi, M., Mateu, J,. and Comas, C. (2020) Directional analysis for point patterns on linear networks. Stat.
#' 
#' @examples
#' # generate random relaisations
#' X <- rpoislpp(0.001,branchnet)
#' M <- angle.lpp(X)


#' @import spatstat
#' @import stats
#' @export
angle.lpp <- function(X,degree=TRUE,directed=TRUE){
  
  n <- npoints(X)
  ang <- matrix(nrow = n,ncol = n)
  for (i in 1:n) {
    
    ps <- psp(rep(X$data$x[i],n),
              rep(X$data$y[i],n),
              X$data$x,
              X$data$y,
              window = X$domain$window)
    ang[,i] <- angles.psp(ps,directed=directed)
  }
  
  if(degree){
    ang[ang>0] <- (ang[ang>0]*360/(2*pi))
    ang[ang<0] <- (ang[ang<0]*360/(2*pi))+360
  }
  
  return(ang)
}