#' Title GOM calcuation with restat in case no-convergence.
#'
#' Function that calculates gom scores. Input is a dataset with variables names, lambdas, and test case values
#' 071717- added cond_func function that calculates conditional survival probabilities and life expectancies
#' @param gom : gom object (pre-exsimated gom object)
#' @param max.iter : maximum interation
#' @param start_vec  : starting value
#' @param tol.gom : convergence criteria
#' @param tol.lik  : bound of likelihood
#'
#' @return a vector of c(updated_gom,iteration,restart.num)
#' @export
#'
#' @examples
calc_gom.restart<- function(gom,
                            max.iter=1000,
                            start_vec=c(0.520308489128769,	0.423437395856997,	0.0000000000000001,	0.0562541150142337),
                            tol.gom=0.000000000001,
                            tol.lik=0.000001){

  updated_gom = calc_gom(gom, max.iter, start_vec=start_vec, tol.gom=tol.gom, tol.lik=tol.lik)[[1]]
  iteration = calc_gom(gom, max.iter, start_vec=start_vec, tol.gom=tol.gom, tol.lik=tol.lik)[[2]]
  restart.num=0
  while(any(updated_gom<0) & restart.num<10){

    updated_gom=round((-1)*updated_gom,10)
    newstart=rep(0,4)
    newstart[which(updated_gom==0)[1]]<-0.8

    if (restart.num==1){updat=rep(0.2/3,3)
    }else{
      updat=rep(0.1,3);updat[restart.num%%3 +1]<-0
    }

    if (any(updated_gom==0)){
      newstart[-which(updated_gom==0)[1]]<-updat
    }else{
      newstart=rep(0.2/3,4);newstart[restart.num%%3 +1]<-0
    }

    updated_gom = calc_gom(gom, max.iter, start_vec=newstart, tol.gom=tol.gom, tol.lik=tol.lik)[[1]]
    iteration = calc_gom(gom, max.iter, start_vec=newstart, tol.gom=tol.gom, tol.lik=tol.lik)[[2]]
    restart.num=restart.num+1
  }
  return(c(updated_gom,iteration,restart.num))
}
