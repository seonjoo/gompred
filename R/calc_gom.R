#' Title gom function
#'
#' @param gom : gom object
#' @param max.iter : max interaction
#' @param start_vec : starting
#' @param tol.gom : convergence
#' @param tol.lik : likelihood
#'
#' @return newlist
#' @export
#'
calc_gom<- function(gom,
                    max.iter=1000,
                    start_vec=c(0.520308489128769,	0.423437395856997,	0.0000000000000001,	0.0562541150142337),
                    tol.gom=0.000000000001,
                    tol.lik=0.000001){
  #Initializing start values
  iteration=1
  #R absolute difference pivot vector
  Rabs_bound_old=1


  temp <-  t(apply(gom[, c('SI', 'SII', 'SIII', 'SIV')], 1, function(x) x* start_vec))
  gom$solution_p <-  apply(temp, 1, sum)

  likelihood<- log(gom$solution_p)

  #validation
  tot_cov<-sum(gom$Recode * gom$Testcase)
  #tot_cov>4

  #Loglikelihood vector
  loglk<- sum(gom$Testcase* log(gom$solution_p))
  lik1<- sum(gom$Recode* gom$SI * gom$Testcase * gom$solution_p^(-1)) - tot_cov
  lik2<- sum(gom$Recode* gom$SII * gom$Testcase * gom$solution_p^(-1)) - tot_cov
  lik3<- sum(gom$Recode* gom$SIII * gom$Testcase * gom$solution_p^(-1)) - tot_cov
  lik4<- sum(gom$Recode* gom$SIV * gom$Testcase * gom$solution_p^(-1)) - tot_cov

  lik<- c(loglk, lik1, lik2, lik3, lik4)

  sum_gradient<- sum(lik[-1])

  pivot<- c(1,1,0,1)
  Rabs_gomvec_old=1

  #Loop until converges
  converge= "Not converged"

  repeat{
    #Adjusted Pivot Vector

    vec<- gom$Recode* gom$Testcase * gom$solution_p^(-2)

    t_1<-apply(gom[, c('SI', 'SII', 'SIII', 'SIV')] , 2, function(x) x*vec)

    hessian<-t(t_1) %*% as.matrix(gom[, c('SI', 'SII', 'SIII', 'SIV')])

    #Checked h1 is created correctly
    h1<- matrix(NA, nrow=nrow(hessian), ncol= ncol(hessian))


    h1[1,1]<-ifelse(pivot[1]==1, ifelse(hessian[1,1]>0.000000000001, 1/ hessian[1,1], 0) , hessian[1,1])

    for (c in 2:4){
      h1[1,c]<-ifelse(pivot[1]==1,h1[1,1]*hessian[1,c],hessian[1,c])
    }

    for (r in 2:4){
      h1[r,1]<- ifelse(pivot[1]==1, -hessian[r,1]*h1[1,1], hessian[r,1])
    }

    for (c in 2:4){
      for (r in 2:4){
        h1[r,c]<- ifelse(pivot[1]==1, hessian[r,c]-hessian[r,1]*h1[1,1]*hessian[c,1], hessian[c,r])
      }}


    #may need to change 1 to R5
    norm_gov_vect<- vector()
    norm_gov_vect[1]<- ifelse(h1[1,1]>0.0000000001, start_vec[1]/1,0.0000000000000001)

    for (c in 2:4){
      norm_gov_vect[c]<- ifelse(h1[1,1]>0.0000000001,start_vec[c]/1,start_vec[c]/(1-start_vec[1]))
    }

    #adjusted pivot vec

    adj_pivot_vec<- ifelse(h1[1,1]>0, pivot[1], 0)
    adj_pivot_vec[2:4]<- pivot[2:4]

    #H2 creation - checked
    h2<- matrix(NA, nrow=nrow(hessian), ncol= ncol(hessian))


    h2[2,2]<-ifelse(adj_pivot_vec[2]==1, ifelse(h1[2,2]>0.000000000001, 1/ h1[2,2], 0) , h1[2,2])

    for (j in c(1,3,4)){
      h2[j,2]<- ifelse(adj_pivot_vec[2]==1, -h1[j,2]*h2[2,2], h1[j,2])
    }

    for (i in  c(1,3,4)){
      h2[2,i]<- ifelse(adj_pivot_vec[2]==1, h2[2,2]*h1[2,i], h1[2,i])
    }

    for (j in c(1,3,4)){
      h2[j,1]<- ifelse(adj_pivot_vec[2]==1, h1[j,1]-h1[j,2]*h2[2,2]*h1[2,1], h1[j,1])
    }

    for (j in c(1,3,4)){
      h2[j,3]<- ifelse(adj_pivot_vec[2]==1, h1[j,3]-h1[j,2]*h2[2,2]*h1[3,2], h1[j,3])
    }

    for (j in c(1,3,4)){
      h2[j,4]<- ifelse(adj_pivot_vec[2]==1, h1[j,4]-h1[j,2]*h2[2,2]*h1[4,2], h1[j,4])
    }

    #Updating pivot vector
    adj_pivot_vec2<- vector()
    adj_pivot_vec2[2]<- ifelse(h2[2,2]>0, adj_pivot_vec[2], 0)
    adj_pivot_vec2[c(1,3,4)]<- adj_pivot_vec[c(1,3,4)]

    #Updating gom vector
    norm_gov_vect2<- vector()
    norm_gov_vect2[1]<- ifelse(h2[2,2]>0.0000000001, norm_gov_vect[1], norm_gov_vect[1]+norm_gov_vect[2])
    norm_gov_vect2[2]<- ifelse(h2[2,2]>0.0000000001, norm_gov_vect[2],0.0000000000000001)

    for (c in 3:4){
      norm_gov_vect2[c]<- norm_gov_vect[c]
    }

    ############
    #H3 and H4 #
    ############

    h3<- matrix(NA, nrow=nrow(hessian), ncol= ncol(hessian))
    h4<- matrix(NA, nrow=nrow(hessian), ncol= ncol(hessian))

    #Creating list of hessian matrices
    hlist<- list(h1,h2, h3, h4)

    #Creating list of gom and adjusted pivot vectors
    adj_pivot_vec3<- vector()
    adj_pivot_vec4<- vector()
    norm_gov_vect3<- vector()
    norm_gov_vect4<- vector()

    pivot_vec_lst<- list(adj_pivot_vec, adj_pivot_vec2, adj_pivot_vec3, adj_pivot_vec4)
    gom_vec_lst<- list(norm_gov_vect, norm_gov_vect2,norm_gov_vect3, norm_gov_vect4)

    colnums<- c(1,2,3,4)

    for (p in c(3:4)){
      hlist[[p]][p,p]<-ifelse(pivot_vec_lst[[p-1]][p]==1, ifelse(hlist[[p-1]][p,p]>0.000000000001, 1/ hlist[[p-1]][p,p], 0) , hlist[[p-1]][p,p])

      remain<- colnums[!colnums %in% p]
      for (j in remain){
        hlist[[p]][j,p]<- ifelse(pivot_vec_lst[[p-1]][p]==1, -hlist[[p-1]][j,p]*hlist[[p]][p,p], hlist[[p-1]][j,p])
      }

      for (i in  remain){
        hlist[[p]][p,i]<- ifelse(pivot_vec_lst[[p-1]][p]==1, hlist[[p]][p,p]*hlist[[p-1]][p,i], hlist[[p-1]][p,i])
      }

      for (l in remain){
        for (j in remain){
          hlist[[p]][j,l]<- ifelse(adj_pivot_vec[p]==1, hlist[[p-1]][j,l]-hlist[[p-1]][j,p]*hlist[[p]][p,p]*hlist[[p-1]][p,l], hlist[[p-1]][j,l])
        }}

      #Computing pivot and updated gom vector for h3
      if (p==3){
        pivot_vec_lst[[p]][p]= ifelse( hlist[[p]][p,p]>0.0000000001,pivot_vec_lst[[p-1]][p],0.0000000000000001)
        for (k in remain){
          pivot_vec_lst[[p]][k]= pivot_vec_lst[[p-1]][k]}

        gom_vec_lst[[p]][p]= ifelse(hlist[[p]][p,p]>0.0000000001, gom_vec_lst[[p-1]][p],0.0000000000000001)

        for (o in 1:2){
          gom_vec_lst[[p]][o]= ifelse(hlist[[p]][p,p]>0.0000000001, gom_vec_lst[[p-1]][o], gom_vec_lst[[p-1]][o]+ hlist[[p-1]][o,p]* gom_vec_lst[[p-1]][p])
        }

        gom_vec_lst[[p]][4]=gom_vec_lst[[p-1]][4]
      }



      #Computing pivot and updated gom vector for h4
      if (p==4){
        pivot_vec_lst[[p]][p]= ifelse( hlist[[p]][p,p]>0,pivot_vec_lst[[p-1]][p],0)
        for (k in remain){
          pivot_vec_lst[[p]][k]= pivot_vec_lst[[p-1]][k]}

        gom_vec_lst[[p]][p]= ifelse(hlist[[p]][p,p]>0.0000000001, gom_vec_lst[[p-1]][p],0.0000000000000001)

        for (o in 1:3){
          gom_vec_lst[[p]][o]= ifelse(hlist[[p]][p,p]>0.0000000001, gom_vec_lst[[p-1]][o], gom_vec_lst[[p-1]][o]+ hlist[[p-1]][o,p]* gom_vec_lst[[p-1]][p])
        }

      }
    }

    #Creating H inverse
    Hinv<-pivot_vec_lst[[4]] %*% t(pivot_vec_lst[[4]]) * hlist[[4]]

    sum_Hinv<- sum(Hinv)
    colsum_Hinv<- apply(Hinv, 2, sum)

    delta<- Hinv %*% lik[2:5]
    C1<- sum(delta)/ sum_Hinv


    epsilon<- colsum_Hinv*C1

    delta_ep<- delta- epsilon

    one_trai_gom<-  gom_vec_lst[[4]] + delta_ep
    adjust_d_e<- ifelse(delta_ep<0, 1-ifelse(one_trai_gom<0, one_trai_gom/ delta_ep, 0), ifelse(one_trai_gom>1, (1-gom_vec_lst[[4]] )/ delta_ep,1))

    R<- min (adjust_d_e)
    two_trai_gom<- gom_vec_lst[[4]] + delta_ep* R

    interior_and_r1<- ifelse(pivot_vec_lst[[4]]==1 & gom_vec_lst[[4]]>0.0001 & gom_vec_lst[[4]]<0.9999 & R>0.0001,1,0)
    trail_gt_small<- ifelse(interior_and_r1==1 & two_trai_gom>=1-0.000000000001,0.99,1)
    trail_lt_small<- ifelse(interior_and_r1==1 & two_trai_gom<=0.000000000001,0.99,1)
    new_mult<- min(trail_gt_small, trail_lt_small)

    Rabs_bound<- R* new_mult

    bound_ind<- ifelse(two_trai_gom>0.00000000000001 & two_trai_gom<0.99999999999999, 0,pivot_vec_lst[[4]])

    updated_gom<- sapply(gom_vec_lst[[4]]+delta_ep*Rabs_bound, function(x) max(x, 0.0000000000000001))
    ll_adj_pivot<- ifelse(updated_gom>0.00000000000001,pivot_vec_lst[[4]],0)

    temp_new <-  t(apply(gom[, c('SI', 'SII', 'SIII', 'SIV')], 1, function(x) x* updated_gom))
    gom$solution_p_new <-  apply(temp_new, 1, sum)



    #gom$solution_p_new <-  t(apply(gom[, c('SI', 'SII', 'SIII', 'SIV')], 1, function(x) x* updated_gom))
    likelihood_new<- log(gom$solution_p_new)

    new_loglk<- sum(gom$Testcase* likelihood_new)

    new_lik1<- sum(gom$Recode* gom$SI * gom$Testcase * gom$solution_p_new^(-1)) - tot_cov
    new_lik2<- sum(gom$Recode* gom$SII * gom$Testcase * gom$solution_p_new^(-1)) - tot_cov
    new_lik3<- sum(gom$Recode* gom$SIII * gom$Testcase * gom$solution_p_new^(-1)) - tot_cov
    new_lik4<- sum(gom$Recode* gom$SIV * gom$Testcase * gom$solution_p_new^(-1)) - tot_cov

    new_lik<- c(new_loglk, new_lik1, new_lik2, new_lik3, new_lik4)

    diff_ll<-new_lik[1] - lik[1]

    #Rabs values at the bottom of sheet
    Rabs_gomvec_new<- sum(abs(updated_gom- start_vec))


    #R abs at top of sheet- uses old numbers
    #Rabs_loglk_grad= sum(new_lik[2:5][new_lik[2:5]>0.000000000001])/ tot_cov
    Rabs_loglk_grad_old= sum(lik[2:5][lik[2:5]>0.000000000001])/ tot_cov
    Rabs_loglk_grad_new= sum(new_lik[2:5][new_lik[2:5]>0.000000000001])/ tot_cov

    #gomcheck1<-ifelse(updated_gom<=0.000000000001 & lik[2:5]>0.000001, 0, 1)
    #gomcheck2<-ifelse(1-updated_gom<=0.000000000001 & lik[2:5]<=-0.000001, 0, 1)

    ###############################
    #print(paste('start_vec','lik','diff_ll','tot_cov'))
    #print(start_vec)
    #print(lik[2:5])
    #print(diff_ll)
    #print(tot_cov)
    ###############################

    gomcheck1<-ifelse(start_vec<=tol.gom & lik[2:5]>tol.lik, 0, 1)
    gomcheck2<-ifelse(1-start_vec<=tol.gom & lik[2:5]<= (-1*tol.lik), 0, 1)

    converge_num<- ifelse(sum(gomcheck1, gomcheck2)<8, 0,1)

    converge<- ifelse(diff_ll/ tot_cov>=-0.00000000000001 & converge_num==1,
                      ifelse( sum(Rabs_gomvec_new, Rabs_loglk_grad_new, diff_ll/tot_cov) <0.000000000001, "Converged","Not converged"), "Not converged")

    if(is.na(converge)){
      updated_gom=c(0.8,0,0,0)
      converge='Not converge'
    }
    # Replacing old information with new

    if (iteration==1) cat("Iteration")
    if (iteration %% 200 ==0)  cat(" ",iteration)

    ####################
    #print(diff_ll/tot_cov)
    ####################
    if (converge=="Converged") {
      cat(": ",converge,"\n")
    }


    last_pivot<- ifelse(updated_gom>0.00000000000001,pivot_vec_lst[[4]],0)

    pivot<- ifelse(last_pivot==1,1,ifelse(new_lik[2:5]<0.000001,0,ifelse(Rabs_gomvec_new>0.1,0,ifelse(Rabs_bound>0.000000000001,1,0))))

    #  if (converge=="Converged") print(converge)


    gom$solution_p= gom$solution_p_new
    start_vec= updated_gom
    lik<- new_lik
    Rabs_gomvec_old= Rabs_gomvec_new


    iteration= iteration +1


    if (iteration==max.iter){
      print("ERROR - COULD NOT CONVERGE")
      updated_gom=updated_gom*(-1)
    }

    if (converge=='Converged' | iteration==max.iter){ break  }
  }

  newlist <- list(updated_gom, iteration)
  return(newlist)

}
