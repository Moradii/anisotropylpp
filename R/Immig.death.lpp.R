#' Simulate immigration-death pattern over  a linear network
#'
#' Simulate immigration-death pattern over  a linear network
#'
#'
#'
#' @param L a linear network
#' @param a anisotropy angle
#' @param t pi/t is the tolerance
#' @param alpha immigration arte
#' @param d death rate
#' @param prob_accep probability to accept a new arriving point
#' @param r r
#' @param N N
#' @param strauss strauss
#' @param burnin burnin
#' @param nsim number of patterns to be simulated
#' @param verbose logical
#' @param Figure logical
#' @param C C
#' 
#' @examples 
#'  X <- Immig.death.lpp(L=branchnet,a=pi/2,t=10,alpha=1.0,d=0.001,prob_accep=0,r=0.15,N=0,strauss=FALSE,burnin=10000,nsim=1,verbose=TRUE,Figure=TRUE,C=1)
#'  plot(X,pch=20,col="red",cols="blue",main="")

#' @references Moradi, M., Mateu, J,. and Comas, C. (2020) Directional analysis for point patterns on linear networks. Stat.
#' @return if \code{nsim=1} a single realization, otherwise a list of realizations of lpp objects
#' @author Carles Comas \email{carles.comas@udl.cat}, Mehdi Moradi \email{m2.moradi@yahoo.com}


#' @import spatstat
#' @import stats
#' @export

Immig.death.lpp <- function(L,a=pi/2,t=10,alpha=1.0,d=0.001,prob_accep=1,r=0.2,N=0,strauss=FALSE,burnin=10000,nsim=1,verbose=TRUE,Figure=TRUE,C=1){
  
  REPS <- burnin
  
  if(nsim>1){
    out <- lapply(X=1:nsim, function(i){
      Function_Immig_death_lpp(L=L,a=a,t=t,alpha=alpha,d=d,prob_accep=prob_accep,r=r,strauss=strauss,burnin=burnin,nsim=1,verbose=TRUE,Figure=TRUE)
    })
  }
  
  
  win <- L$window
  
    Ys <- rpoislpp(100/volume(L),L) #define a newly arrived point
   Xap <- as.numeric(Ys$data$x[1])
   Yap <- as.numeric(Ys$data$y[1])
  pppi <- lpp(ppp(as.numeric(Xap),as.numeric(Yap),window=win),L)
  
  nr <- npoints(pppi)
  # death<-c()
  death <- rep(0,nr)
  nt <- c()
  
  int <- alpha/d*C
  if(Figure){
    plot(1,nr,xlim=c(0,REPS),ylim=c(0,int)) 
  }
  
  
  for (i1 in 1:REPS){
    nr <- npoints(pppi)
    nt[i1] <- nr
    
    if(N>0){
      if(i1 > round(0.9*REPS) & nr==N){break}
    }
    
    if (verbose) {
      cat(paste(i1), ",") #cat(paste(i1), "\n")
      flush.console()
    }
    
    
    prob_set <- 1 #iNITIALLY THIS PROBABLITY IS EQUAL TO 1
    
    #Test for death of a point
    
    if(nr>0){
      for(i in 1:nr){
        Z <- runif(1,0,1)
        if(Z<=d){ #delete a point at random
          death[i] <- 1
        }
      }
      
      j=0
      xdp <- c(); ydp <- c(); death2 <- c()
      for(i in 1:nr){
        if(death[i]==0){
          j <- j+1
          xdp[j] <- pppi$data$x[i]
          ydp[j] <- pppi$data$y[i]
          death2[j] <- death[i]
        }
      }
      death  <-death2
      pppi <- lpp(ppp(as.numeric(xdp),as.numeric(ydp),window=win),L)
      nr <- length(pppi$data$x)
    } #if(nr>0){
    #TEST FOR INMIGRANS
    #HERE IS TO DECIDE TO ACCEPT O REJECT A NEWLY ARRIVED INMIGRANT
    Z<-runif(1,0,1)
    
    if(Z<=alpha){      
      prob_set <- 1 #INITIALLY THIS PROBABLITY IS EQUAL TO 1
            Ys <- rpoislpp(100/volume(L),L) #define newly arrived point
           Xap <- as.numeric(Ys$data$x[1])
           Yap <- as.numeric(Ys$data$y[1])
      
            ps <- psp(as.numeric(rep(Xap,nr)),as.numeric(rep(Yap,nr)),as.numeric(pppi$data$x),as.numeric(pppi$data$y),
                window = pppi$domain$window,check=TRUE)
           ang <- angles.psp(ps)
         pppi2 <-pppi
      
        #pppi2 <- pppi2[ang <(a+pi/t) & ang > (a-(pi/t))]
      
         pppi2 <- pppi2[ang < (a-(pi/t)) | ang >(a+pi/t)]
      
            Ys <- lpp(ppp(Xap,Yap,window=win),L)
      
            Xs <- superimpose(pppi2,Ys) #new lpp
            nr <- npoints(pppi2)
      
      if(strauss){
        Xs <- superimpose(pppi,Ys)
        nr <- npoints(pppi)
      }
      
           dist <- array(0,c(nr+1,nr+1))
      dim(dist) <- c(nr+1,nr+1)
           dist <- pairdist(Xs) ##the vector with the pairwise distance between points of Xs1
      if(nr>0){
        for(i in 1:nr){
          DIST_s <- dist[i,nr+1]
          if(DIST_s <= r){
            #COMPUTED THE SPATIAL PROJECTED ANGLE BETWEEN POINTS (Xap,Yap) i (X[i],Y[i])
            prob_set <- prob_set*prob_accep 
          } #DIST_s<r 
        } #for(i in 1:nr)
      } #if nr>0
      
      #ONCE prob_set IS COMPUTED THEN WE HAVE TO TEST IF THE POINT IS ACCEPTED OR NOT
       Z <- runif(1,0,1)
      Xs <- superimpose(pppi,Ys)
      nr <-npoints(pppi)
      if(Z<=prob_set){
               pppi <- Xs
        death[nr+1] <- 0
      }              
    }#z<=alpha
    if(Figure){
      points(i1,nr)
    }
    
  } #REPS
  
  
  if(nsim==1){
    out <- pppi
     nt <- nt
  }
  
  attr(out,"nt") <- nt
  
  # result<-list(out=out,nt=nt)
  # invisible(return(result))
  
  return(out)
}




