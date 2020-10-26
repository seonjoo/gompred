#' Title Computing TLE
#'
#' @param GOM1
#' @param GOM2
#' @param GOM3
#' @param GOM4
#' @param tle_lam
#'
#' @return list(tle, S_tle)
#' @export
#'
tle<- function(GOM1, GOM2, GOM3, GOM4, tle_lam){

  S_tle<- matrix(NA, nrow=21, ncol=1)
  S_tle[1]=1

  for (j in 2:21){
    S_tle[j]= S_tle[j-1]* ( GOM1*tle_lam[1, j] + GOM2*tle_lam[2, j] + GOM3*tle_lam[3, j] + GOM4*tle_lam[4, j] )
  }


  first<- sum(S_tle[1:20])
  second<- sum(S_tle[2:21])
  tle=round(sum(sum(first), sum(second))/4, digits=2)

  #tle<-print(tle)
  S_tle<- S_tle

  return(list(tle, S_tle))

}
