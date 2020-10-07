#' Simulate anisotropic Matern cluster process on linear networks
#'
#' Simulate anisotropic Matern cluster process on linear networks based on given interaction distance and angle 
#'
#'
#'
#' @param pint intensity function of the parents if model is Poisson, otherwise the number of parrents
#' @param r iteraction distance 
#' @param cint intensity function of offsprings if model is Poisson, otherwise the number of offsprings per cluster
#' @param L a linear network
#' @param a angle of interaction
#' @param t pi/t will be considered as a tolerance for angle of interaction
#' @param model model to be considered for simulating parents and offsprings, currently only Poisson and Uniform. Default is Poisson.
#' @param relaxed logical. if TRUE the function keeps some points which do not fall within angle condition according to a given probability
#' @param prob if relaxed=TRUE, then the simulation also includes points, which do not fall within an angle a with respect to their parent, with probability prob
#' @param nsim number of simulations
#' @param check logical
#' 
#' @details this function generates realizations from (an)isotropic Matern Cluster point processes over linear networks. the strength of anisotropy can be controlled by \code{relaxed} and \code{prob}
#' 
#' @references Moradi, M., Mateu, J,. and Comas, C. (2020) Directional analysis for point patterns on linear networks. Stat.
#' @return if \code{nsim=1} a single realization, otherwise a list of realizations of lpp objects
#' @seealso \link[spatstat]{rMatClust}
#' @author Mehdi Moradi \email{m2.moradi@yahoo.com}
#' @examples
#' # generate random realizations
#' X <- r.anis.MatClust.lpp(pint=0.01,r=300,cint=0.2,L=branchnet,t=50)
#' plot(X,main="",cols=2,col=4,pch=20)
#' Y <- r.anis.MatClust.lpp(pint=0.01,r=300,cint=0.2,L=branchnet,t=50,a=pi/2)
#' plot(Y,main="",cols=2,col=4,pch=20)



#' @import spatstat
#' @import stats
#' @export
r.anis.MatClust.lpp <- function(pint=10,r=0.2,cint=10,L=L,a=pi/3,t=10, model=c("Poisson","Uniform"),
                                relaxed=FALSE,prob=0.5,nsim=1,check=FALSE){
  
  if(nsim>1){
    out <- lapply(X=1:nsim, function(i){
      r.anis.MatClust.lpp(pint = pint,r=r,cint = cint,L=L,a=a,t=t,nsim=1,model = model,
                          relaxed=relaxed,prob=prob,check=check)
    })
    return(out)
  }
  
  xFinal <- rep(0,0)
  yFinal <- rep(0,0)
  
  if(missing(model)) model="Poisson"
  
  if(model=="Poisson"){
    X <- rpoislpp(pint,L)
    
    cl <- vector("list",npoints(X))
    
    h=1
    
    for (i in 1:npoints(X)){
      
      sub <- lineardisc(L,c(X$data$x[i],X$data$y[i]),r,plotit = F)
      sublin <- sub$lines
      sublin <- as.linnet(sublin)
      u <- rpoislpp(cint,sublin)
      
      ps <- psp(rep(X$data$x[i],npoints(u)),
                rep(X$data$y[i],npoints(u)),
                u$data$x,
                u$data$y,
                window = X$domain$window,check=check)
      
      ang <- angles.psp(ps)
      uout <- u[ang <= (a+(pi/t)) & ang >= (a-(pi/t))]
      
      if(relaxed){
        ustor <- u[ang > (a+(pi/t)) | ang < (a-(pi/t))]
        ustorthin <- rthin(ustor,P=prob)
        cl[[i]] <- superimpose(uout,ustorthin)
      }else{
        cl[[i]] <- uout
      }
      
    }
  }else if(model=="Uniform"){
    
    X <- runiflpp(pint,L)
    
    cl <- vector("list",npoints(X))
    
    h=1
    
    for (i in 1:npoints(X)){
      
      sub <- lineardisc(L,c(X$data$x[i],X$data$y[i]),r,plotit = F)
      sublin <- sub$lines
      sublin <- as.linnet(sublin)
      u <- runiflpp(cint,sublin)
      
      ps <- psp(rep(X$data$x[i],npoints(u)),
                rep(X$data$y[i],npoints(u)),
                u$data$x,
                u$data$y,
                window = X$domain$window,check=check)
      
      ang <- angles.psp(ps)
      uout <- u[ang <= (a+(pi/t)) & ang >= (a-(pi/t))]
      
      if(relaxed){
        ustor <- u[ang > (a+(pi/t)) | ang < (a-(pi/t))]
        ustorthin <- rthin(ustor,P=prob)
        cl[[i]] <- superimpose(uout,ustorthin)
      }else{
        cl[[i]] <- uout
      }
      
    }
  }else {stop("model should be either Poisson or Uniform")}
  
  for (j in 1:npoints(X)){
    if(npoints(cl[[j]])>0){
      for(k in 1:length(cl[[j]]$data$x)){
        xFinal[h] <- cl[[j]]$data$x[k]
        yFinal[h] <- cl[[j]]$data$y[k]
        h=h+1
      }
    }
  }
  
  out <- lpp(data.frame(unlist(xFinal),unlist(yFinal)),L)
  
  return(out)
}
