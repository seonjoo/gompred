#' Title  Computing DFLE and DLE
#'
#' @param GOM1
#' @param GOM2
#' @param GOM3
#' @param GOM4
#' @param tle_lam
#' @param dle_lam
#'
#' @return list(results, S_dle, S_nle)
#' @export
#'
#' @examples
dle_nle<- function(GOM1, GOM2, GOM3, GOM4, tle_lam, dle_lam){

  first_part<- matrix(NA, nrow=1, ncol=22)
  second_part<- matrix(NA, nrow=1, ncol=22)
  S_dle<- matrix(NA, ncol=1, nrow=21)
  S_nle<- matrix(NA, ncol=1, nrow=21)

  first_part[1]=1

  for (j in 2:22){
    i= j-1
    first_part[j]= (first_part[j-1]* (GOM1*tle_lam[1, j] + GOM2*tle_lam[2, j] + GOM3*tle_lam[3, j] + GOM4*tle_lam[4, j]))
    second_part[i]= GOM1*dle_lam[1, j] + GOM2*dle_lam[2, j] + GOM3*dle_lam[3, j] + GOM4*dle_lam[4, j]
  }

  for (k in 1:21){
    S_dle[k]= first_part[k] *second_part[k]
    S_nle[k]= first_part[k]- S_dle[k]
  }

  first= sum(S_dle[1:20])
  second=sum(S_dle[2:21])


  dle=sum(sum(first), sum(second))/4
  nle=sum(sum(S_nle[1:20]), sum(S_nle[2:21]))/4

  results= cbind(round(dle, digits=2), round(nle, digits=2))

  return(list(results, S_dle, S_nle))
}
