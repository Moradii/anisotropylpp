#' Inhomogeneous K-function accounting for directions
#'
#' Inhomogeneous K-function accounting for directions
#'
#'
#'
#' @param X an object of class lpp
#' @param lambda estimated intensity at data points. If not given, it will be automatically obtained from \link[spatstat]{densityQuick.lpp}
#' @param normalize logical. whether to normalize the estimate. See \link[spatstat]{linearKinhom}
#' @param r optional. distance vector wherein the K-funtion is evaluated
#' @param phi optional. angle vector wherein the K-funtion is evaluated
#' @param maxphi max angle. this will be used only when phi is not given
#' @param edge optional. whether to use edge correction or not
#' @param nxy dimension. this will be considered as the length of r and phi when they are not given
#' @param parallel logical. if TRUE, it uses \link[parallel]{mclapply}.
#' @param verbose logical. if TRUE, it prints the progress of the function.
#' @param ... arguments passed to \link[parallel]{mclapply} for configuration of the cores.
#' 
#' @details This function calculates the inhomogeneous K-function over a grid of distances and angles.
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' @seealso \link[spatstat]{linearKinhom}
#' @return an object of class sumlpp
#' @note this function can be quite slow, we suggest to use \link[anisotropylpp]{nnangle.lpp}
#' @references Moradi, M., Mateu, J,. and Comas, C. (2020) Directional analysis for point patterns on linear networks. Stat.
#' @examples
#' # generate random relaisations
#' X <- rpoislpp(0.001,branchnet)
#' K <- anisLinearKinhom(X)


#' @import spatstat
#' @import stats
#' @export

anisLinearKinhom <- function(X,lambda=lambda,normalize=TRUE,r=NULL,
                             phi=NULL,maxphi= 360,edge=TRUE,nxy=10,parallel=FALSE,
                             verbose=FALSE,...){
  if (!inherits(X, "lpp")) stop("X should be from class lpp")
  
  if(missing(lambda)) {
    
    lambda <- densityQuick.lpp(X,sigma=bw.scott.iso(X),at="points")
  }
  
  if(length(lambda)!=npoints(X)) stop("lambda should be a vector of intensity values at data points")
  
  l <- domain(X)
  sdist <- pairdist.lpp(X)
  phiang <- angle.lpp(X)
  
  n <- npoints(X)
  
  lamden <- outer(lambda,lambda,FUN = "*")
  diag(lamden) <- 1
  
  maxs <- 0.7*max(sdist[!is.infinite(sdist)])
  
  if(is.null(r) & is.null(phi)){
    phi <- seq(0,maxphi,length.out = nxy)
    r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
    K <- matrix(NA,nrow = nxy,ncol = nxy)
  }else{
    if(is.null(r)){
      r <- seq((maxs/nxy),maxs,by=(maxs-(maxs/nxy))/(nxy-1))
      K <- matrix(NA,nrow = length(r),ncol = length(phi))
    }
    else if(is.null(phi)){
      phi <- seq((maxphi/nxy),maxphi,by=(maxphi-(maxphi/nxy))/(nxy-1))
      K <- matrix(NA,nrow = length(r),ncol = length(phi))
    }else{
      K <- matrix(NA,nrow = length(r),ncol = length(phi))
    }
    
  }
  
  

  if(edge){
    
    
    if(parallel){
      
      edge1 <- mclapply(X=1:n, function(j){
        
        if(verbose){
          if(j<n) {
            cat(paste(j),",")
            flush.console()
          } else {
            cat(paste(j),"\n")
            flush.console()
          }
        }
        
        lapply(X=1:(n-1), function(i){
          ldisc <- lineardisc(l, X[-j][i], sdist[-j,j][i], plotit=FALSE)
          sup <- superimpose(X[-j][i],ldisc$endpoints)
          ang <- angle.lpp(sup)[-1,1]
          list(ang,(npoints(ldisc$endpoints)^2))
        })
        
      },...)
      
      
    }else{
      
      edge1 <- lapply(X=1:n, function(j){
        
        if(verbose){
          if(j<n) {
            cat(paste(j),",")
            flush.console()
          } else {
            cat(paste(j),"\n")
            flush.console()
          }
        }
        
        lapply(X=1:(n-1), function(i){
          ldisc <- lineardisc(l, X[-j][i], sdist[-j,j][i], plotit=FALSE)
          sup <- superimpose(X[-j][i],ldisc$endpoints)
          ang <- angle.lpp(sup)[-1,1]
          list(ang,(npoints(ldisc$endpoints)^2))
        })
        
      })
      
    }
    
  }
  
  
  
  
  for (h in 1:length(r)){
    
    for (v in 1:length(phi)) {
      
      if(edge){
        
        ml <- matrix(1, n, n)
        
        for (i in 1:n) {
          edgein <- c()
          for (j in 1:(n-1)) {
            edgein[j] <- sum(edge1[[i]][[j]][[1]]<=phi[v])/edge1[[i]][[j]][[2]]
          }
          ml[ -i, i] <- edgein
        }
      }
      
      out <- (sdist <= r[h])*(phiang <= phi[v])
      diag(out) <- 0
      
      if(edge){
      
        kout <- out*ml/lamden
        
      }else{
        
        kout <- out/lamden
        
      }
      
      K[h,v] <- sum(kout[!is.na(kout) & !is.infinite(kout)])
    }
  }
  
  if(normalize){
    revrho <- outer(1/(lambda),1/(lambda),FUN = "*")
    appx <- (volume(l))/(sum(revrho[lower.tri(revrho, diag = FALSE)])*2)
    K <- K*appx
  }
  else{
    K <- K/(volume(l))
  }
  
  
  Kout <- list(Kest=K,r=r,phi=phi)
  
  class(Kout) <- c("sumlpp")
  return(Kout)
}

print.sumlpp <- function(x){
  cat(paste0("inhomogeneous K-function"),"\n")
}
