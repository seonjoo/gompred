#' Title Compute conditional survival functions
#'
#' @param GOM1
#' @param GOM2
#' @param GOM3
#' @param GOM4
#' @param tle_lam
#' @param dle_lam
#'
#' @return cond_func : returns
#' @export
#'
cond_func<- function(GOM1, GOM2, GOM3, GOM4, tle_lam, dle_lam){

  #Bomb if GOM scores don't sum to 1 - validation done in shiny app
  #stopifnot(sum(GOM1, GOM2, GOM3, GOM4)==1)

  S_tle<- matrix(NA, nrow=21, ncol=1)
  pi<- matrix(NA, nrow=21, ncol=1)
  ppi<-matrix(NA, nrow=21, ncol=1)
  ppic<-matrix(NA, nrow=22, ncol=1)
  dd<- matrix(NA, nrow=21, ncol=1)
  hp<- matrix(NA, nrow=21, ncol=1)
  d<- matrix(NA, nrow=21, ncol=1)
  f<- matrix(NA, nrow=21, ncol=1)

  S_tle[1]=1
  ppi[1]=1
  ppic[1]=1

  for (j in 2:22){
    S_tle[j]= S_tle[j-1]* ( GOM1*tle_lam[1, j] + GOM2*tle_lam[2, j] + GOM3*tle_lam[3, j] + GOM4*tle_lam[4, j] )
  }

  for (j in 1:21){
    pi[j]<- GOM1*dle_lam[1, j+1] + GOM2*dle_lam[2, j+1] + GOM3*dle_lam[3, j+1] + GOM4*dle_lam[4, j+1]
  }

  for (i in 2:22){
    ppi[i]<- S_tle[i-1]*pi[i-1]
  }

  for (k in 1:21){
    ppic[k+1]<- S_tle[k]-ppi[k+1]
  }

  for (m in 1:21){
    dd[m]= ppic[m]- ppic[m+1]
  }

  ppic[23]= ppic[22]^2/ppic[21]

  for (n in 1:21){
    hp[n]= log(S_tle[n+1]/ S_tle[n])/log(ppic[n+2]/ppic[n+1])
  }

  d[1]= dd[1]
  for (n in 2:21){
    d[n]= dd[n]*(1-hp[n-1])
  }
  l20= ppic[20]
  d[22]= l20*(1- hp[21])

  sum_d= sum(d[1:21])


  f[1]= (d[1]+ d[1+1]/2)/sum_d
  for (r in 2:21){
    f[r]= (d[r]/2+ d[r+1]/2)/sum_d
  }

  IS= 1/S_tle
  IS=IS[1:21]

  #Survival with FTC
  S_ftc= matrix(NA, nrow=21, ncol=1)
  S_ftc[1]=1

  for (b in 2:21){
    S_ftc[b]=sum(f[1:(20-b+2)]* IS[1:(20-b+2)] * S_tle[b:21])
  }

  terminal_adj= 0.5*(S_ftc[1] + (sum(f[1:20]*IS[1:20])*S_tle[21]) )

  #Computing normalized survival without FTC
  S_wo_ftc<- matrix(NA, nrow=21, ncol=1)
  S_wo_ftc[1]=1
  for (b in 2:21){
    S_wo_ftc[b]= ppic[b+1]/ ppic[2]
  }

  #Computing decrements from normalized survival without FTC
  decr_wo_ftc<- matrix(NA, nrow=21, ncol=1)
  for (b in 2:21){
    decr_wo_ftc[b-1]= S_wo_ftc[b-1] - S_wo_ftc[b]
  }
  decr_wo_ftc[21]=  S_wo_ftc[21]

  #Sum of decrements w/o FTC
  sum_decr_wo_ftc= sum(decr_wo_ftc)

  #Death decrements wo FTC
  decr_death_wo_ftc<- matrix(NA, nrow=21, ncol=1)
  for (b in 1:21){
    decr_death_wo_ftc[b]= decr_wo_ftc[b]*hp[b]
  }

  #Computing time sequence
  ts= seq(from=0, to=20, by=1)

  #Computing hazard proportion ts
  hazard_ts= unlist(lapply(ts, function(x) x/2 +0.25))
  hazard_ts[21]= ts[21]/2

  #Computing FTC decrements from normalized surivval without FTC
  decr_ftc_wo_ftc<- decr_wo_ftc-decr_death_wo_ftc

  #Computing SDLE
  Sdle= (sum(S_ftc)- terminal_adj)/2

  #Computing cdle_2
  #Have to find tle
  tle=tle(GOM1, GOM2, GOM3, GOM4, tle_lam)[[1]]
  pi_0_tle=pi[1]*tle
  p_ftc_n= sum(d)-pi[1]

  cdle_2=(Sdle-pi_0_tle/sum(d))*sum(d)/p_ftc_n
  #Computing Time-to-FTC
  ttftc<-sum(decr_ftc_wo_ftc*hazard_ts)/sum(decr_ftc_wo_ftc)


  #Computing conditional dlfe
  #Have to get NLE from dle_nle function
  nle=dle_nle(GOM1, GOM2, GOM3, GOM4, tle_lam, dle_lam)[[1]][2]
  cond_nle= nle/ ppic[2]

  #Computing Conditional NLE variable for input to time-to-death variable
  cond_nle_temp=  sum(hazard_ts*decr_death_wo_ftc)

  #Computing time to death conditional on no ftc at intake
  ttd_noftc= cond_nle_temp/sum(decr_death_wo_ftc)

  #Computing total time to death conditional on no ftc at intake
  tot_time_death<- cdle_2 + ttftc

  #Creating a vector of conditional life expectancies
  cond_exp<- c(cond_nle, ttd_noftc, ttftc, tot_time_death,  Sdle, tle, cdle_2)
  names(cond_exp)<- c('exp2', 'exp3', 'exp4', 'exp5', 'exp6', 'exp7', 'exp8')

  #Creating vector of conditional probabilities
  #Computing total prob- no FTC at intake
  prob2<- 1- ppi[2]
  #Computing prob death before FTC
  prob3<- prob2- p_ftc_n
  #Computing probability of FTC before death
  prob4<-  p_ftc_n
  #Computing prob of death after ftc
  prob5<-  p_ftc_n
  prob6<- sum_d
  prob7<-ppi[2]
  prob8<- p_ftc_n

  cond_prob<- c(prob2, prob3, prob4, prob5, prob6, prob7, prob8)
  names(cond_prob) <- c('prob2', 'prob3', 'prob4', 'prob5', 'prob6', 'prob7', 'prob8')

  #Outputting list with conditional probabilites and life expectancies
  final_cond<- list(cond_exp, cond_prob*100)

  #Rounding output
  final_cond[[1]]<-round(final_cond[[1]], 2)
  final_cond[[2]]<-round(final_cond[[2]], 1)

  # output survival for FTC and normalized survival
  final_cond[[3]]<-S_ftc
  final_cond[[4]]<-S_wo_ftc
  final_cond[[5]]<-pi

  return(final_cond)
}
