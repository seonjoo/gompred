#' Title Com computation based on the input data and gom object.
#'
#' @param inputdata
#' @param inputdata_na
#' @param gom.obj
#' @param out.excel.file
#' @param start_vec1
#' @param start_vec2
#' @param max.iter
#' @param tol.gom
#' @param tol.lik
#'
#' @import shiny
#' @import ggplot2
#' @import openxlsx
#' @import tidyverse
#' @import readxl
#' @export
#'
gom.qing<-function(inputdata,inputdata_na, #main dataset and NA dataset that need a diff start vector
                   gom.obj, ## gom object that Eric provides
                   out.excel.file, ## excel file name
                   start_vec1,start_vec2, # two differenct start vectors for data1 and data2
                   max.iter=1000, ## maximum iteration
                   tol.gom=10^(-12),
                   tol.lik=10^(-6)

){

  ################
  #SET UP  ###
  ################
  testdata<- inputdata

  #create excel workbook
  GOM_results <- createWorkbook()
  addWorksheet(GOM_results, sheetName = "GOM Scores")
  addWorksheet(GOM_results, sheetName = "Life Expectancy")
  addWorksheet(GOM_results, sheetName = "Survival")
  addWorksheet(GOM_results, sheetName = "Percents Table")
  addWorksheet(GOM_results, sheetName = "Cond1 Table")
  addWorksheet(GOM_results, sheetName = "Cond2 Table")

  #create data frames
  gom_scores_all = data.frame()
  le_all = data.frame()
  survival_all = data.frame()
  percents_all = data.frame()
  cond1_all = data.frame()
  cond2_all = data.frame()


  ################
  #GOM scores  ###
  ################


  for (i in 1:nrow(testdata)){ ### SL Comment 1. Change to nrow(testdata)
    testcase <- testdata[i,]
    id <- testdata[i,1]
    print(id)

    if (id %in% inputdata_na$id){
      start_vec = start_vec2
    }else{
      start_vec = start_vec1
    }
    ### This section reoder variables to make the input files comparabiel to the gom function
    ### missing values will have NA value
    selection <- lapply(choice, function(x) x[[1]])
    selection[names(testcase)] <- testcase
    selection

    #Creating dataset from input values
    out<- as.data.frame(data.table::transpose(selection))
    rownames(out) <- colnames(as.data.frame(selection))
    #colnames(out)<- 'Resp Number'
    colnames(out)<- 'Resp'
    out$`Var Name`<- rownames(out)
    out$`Testcase`<- 1

    #Merging gom data and patient data
    gom.obj$orderg<- 1:nrow(gom.obj)
    #final<- merge(gom.obj, out, by=c("Var Name", "Resp Number"), all.x=TRUE)
    final<- merge(gom.obj, out, by=c("Var Name", "Resp"), all.x=TRUE)
    final<- final[order(gom.obj$orderg), ]
    final$Testcase[is.na(final$Testcase)]<-0
    #sum(final$Testcase)

    validate(
      #fix spelling error provied s/b provided -- Eric 11/06/18
      need( sum(final$Testcase)> 4, "Not enough information provided to compute Gom Scores")
    )

    tmp <-NA

    if ( sum(final$Testcase)> 4 ){
      try(
        tmp <- t(round(calc_gom.restart(final,max.iter=max.iter,start_vec=start_vec,tol.gom=tol.gom,tol.lik=tol.lik), digits =9))
      )
    }

    # if gom_score = na then also output
    #while (all(is.na(tmp)) | any(tmp<0))
    #{
    #  tmp <- t(round(calc_gom.restart(final,max.iter=max.iter,start_vec=c(rep(4,0,1)),tol.gom=tol.gom,tol.lik=tol.lik), digits =9))
    #}


    if (all(is.na(tmp)) | any(tmp<0)){
      gom_scores = t(c(888,888,888,888))
      iteration = tmp[5]
      restart.num=NA
    }else{
      gom_scores=t(tmp[1:4])
      iteration = tmp[5]
      restart.num=tmp[6]
    }

    #add id
    gom_scores_id<-data.frame(ID=id,GOM=gom_scores,Iteration = iteration, Num_Cov=sum(final$Testcase),restart.num=restart.num)

    #combine gom_scores for all id
    gom_scores_all<-  rbind(gom_scores_all, gom_scores_id)

    #######################################
    #GOM scores to compute TLE, DFLE, DLE, FTC #
    #######################################
    #Creating plot dataset

    dataset<-gom_scores
    time<-seq(0, 10, 0.5)
    S_tle<-unlist(tle(dataset[1], dataset[2], dataset[3], dataset[4], tle_lam)[[2]])
    S_dle<- dle_nle(dataset[1], dataset[2], dataset[3], dataset[4],  tle_lam, dle_lam)[[2]]
    S_nle<- dle_nle(dataset[1], dataset[2], dataset[3], dataset[4],  tle_lam, dle_lam)[[3]]

    S_ftc<-cond_func(gom_scores[1], gom_scores[2], gom_scores[3], gom_scores[4], tle_lam, dle_lam)[[3]]
    S_wo_ftc<-cond_func(gom_scores[1], gom_scores[2], gom_scores[3], gom_scores[4], tle_lam, dle_lam)[[4]]
    pi<-cond_func(gom_scores[1], gom_scores[2], gom_scores[3], gom_scores[4], tle_lam, dle_lam)[[5]]

    survival_tle<-cbind(time, as.data.frame(S_tle), 'TLE')
    colnames(survival_tle)[2:3]<- c('S', 'Outcome')

    survival_dle<- cbind(time, as.data.frame(S_dle), 'DLE')
    colnames(survival_dle)[2:3]<- c('S', 'Outcome')


    survival_nle<- cbind(time, as.data.frame(S_nle), 'DFLE')
    colnames(survival_nle)[2:3]<- c('S', 'Outcome')

    survival_ftc<-cbind(time, as.data.frame(S_ftc), 'FTC')
    colnames(survival_ftc)[2:3]<- c('S', 'Outcome')

    survival_wo_ftc<-cbind(time, as.data.frame(S_wo_ftc), 'Wo_FTC')
    colnames(survival_wo_ftc)[2:3]<- c('S', 'Outcome')

    survival_pi<-cbind(time, as.data.frame(pi), 'pi')
    colnames(survival_pi)[2:3]<- c('S', 'Outcome')
    #Combining


    survival_plot<- rbind(survival_tle, survival_dle, survival_nle,survival_ftc,survival_wo_ftc,survival_pi)
    # transpose
    output_survival_plot<- spread(survival_plot,time, S)
    # add id
    survival_id<-cbind(id,output_survival_plot)

    #combine gom_scores for all id
    survival_all<-  rbind(survival_all, survival_id)

    #Output
    output_tle <- tle(dataset[1], dataset[2], dataset[3], dataset[4], tle_lam)[[1]]
    output_dle <- dle_nle(dataset[1], dataset[2], dataset[3], dataset[4], tle_lam, dle_lam)[[1]][1]
    output_nle <-dle_nle(dataset[1], dataset[2], dataset[3], dataset[4], tle_lam, dle_lam)[[1]][2]

    output_le<- cbind(output_tle,output_dle, output_nle)

    output_le_id<-cbind(id,output_le)

    #combine le for all id
    le_all<-  rbind(le_all, output_le_id)



    #######################################
    #GOM scores to compute survival plot #
    #######################################


    #output_plot<-ggplot( survival_plot , aes(x=time, y=S, color= Outcome, shape=Outcome)) +
    #  geom_point() + geom_line() + ggtitle('') +
    #  xlab('Years since Intake') + ylab('Probability of Survival') +
    #  scale_x_continuous(breaks=seq(0,10,1), minor_breaks=NULL, expand=c(0,0)) +
    #  scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=NULL, expand=c(0,0)) +
    #  theme(text = element_text(size=20)) +
    #  ggtitle(paste("ID = ", id))
    #print(output_plot)

    # insertPlot(GOM_results, sheet = "PLOTs",  xy=c("I", 20*i-19))

    ################
    # Percents   ###
    ################

    #Getting corresponding values using linear interpolation
    f <- approxfun(survival_tle$S, survival_tle$time)
    corr_values= as.data.frame(matrix(NA, ncol=2, nrow=5))
    colnames(corr_values)= c('Survival Quantile', 'Years')
    corr_values[1,]= c(90, f(.9))
    corr_values[2,]= c(75, f(.75))
    corr_values[3,]= c(50, f(.5))
    corr_values[4,]= c(25, f(.25))
    corr_values[5,]= c(10, f(.1))

    #Finding the last percentile in 10 years
    minpercent<-min(corr_values$`Survival Quantile`[!is.na(corr_values$Years)])

    if (minpercent==10) {
      final_prob_list= corr_values
    }

    f2 <- approxfun(survival_tle$time, survival_tle$S)
    if (minpercent!=10) {
      Y10= c(f2(10)*100, 10)
      final_prob_list= rbind(corr_values[!is.na(corr_values$Years),], Y10)
    }
    #increase digits to 2 -- Eric 11/06/18

    output_percents<- round(final_prob_list[,c(2,1)], digits=2)

    # add id
    output_percents_id<-cbind(id,output_percents)

    #combine percent for all id
    percents_all<-  rbind(percents_all, output_percents_id)

    #####################
    ##  Cond 1 Table  ###
    #####################
    #caption= "10- Year Conditional Life Expectancies"
    cond_exp1<-data.frame(matrix(cond_func(gom_scores[1], gom_scores[2], gom_scores[3], gom_scores[4], tle_lam, dle_lam)[[1]], nrow=1))
    colnames(cond_exp1)<- c('No FTC at Intake: Total (2)', 'No FTC at Intake: Time to Death (3)', 'No FTC at Intake: Time to FTC (4)',
                            'No FTC at Intake: Total Time to Death (5)', 'FTC At any Time: Total (6)', 'FTC at any Time: FTC at Intake (7)',
                            'FTC at Any Time: No FTC Intake; FTC before death (8)')
    output_cond1<-cond_exp1
    # add id
    output_cond1_id<-cbind(id,output_cond1)

    #combine percent for all id
    cond1_all<-  rbind(cond1_all, output_cond1_id)

    #####################
    ##  Cond 2 Table  ###
    #####################
    #caption= "Probabilites of Disability and Death over 10 Years"

    cond_exp2<-data.frame(matrix(cond_func(gom_scores[1], gom_scores[2], gom_scores[3], gom_scores[4], tle_lam, dle_lam)[[2]], nrow=1))
    colnames(cond_exp2)<- c('No FTC at Intake: Total (2)', 'No FTC at Intake: Death before FTC (3)', 'No FTC at Intake: FTC before Death (4)',
                            'No FTC at Intake: Death after FTC (5)', 'FTC At any Time: Total (6)', 'FTC at any Time: FTC at Intake (7)',
                            'FTC at Any Time: No FTC Intake; FTC before death (8)')
    output_cond2<-cond_exp2
    # add id
    output_cond2_id<-cbind(id,output_cond2)

    #combine percent for all id
    cond2_all<-  rbind(cond2_all, output_cond2_id)


  }  # end of for loop


  #colnames(gom_scores_all)<- c("ID","GOM1", "GOM2", "GOM3", "GOM4")
  colnames(le_all)<- c("ID","TLE", "DLE", "NLE")


  ###################
  # output to EXCEL #
  ###################

  #gom scores
  writeDataTable(GOM_results, sheet = "GOM Scores", x = gom_scores_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  #tle, dle, nle
  writeDataTable(GOM_results, sheet = "Life Expectancy", x = le_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  #survival data
  survival_mean_tle<- survival_all%>%
    filter(Outcome == 'TLE')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_mean_dle<- survival_all%>%
    filter(Outcome == 'DLE')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_mean_dfle<- survival_all%>%
    filter(Outcome == 'DFLE')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_mean_ftc<- survival_all%>%
    filter(Outcome == 'FTC')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_mean_wo_ftc<- survival_all%>%
    filter(Outcome == 'Wo_FTC')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_pi<- survival_all%>%
    filter(Outcome == 'pi')%>%
    select(-id,-Outcome)%>%
    colMeans(na.rm = FALSE, dims = 1)

  survival_mean<- as.data.frame(rbind(survival_mean_tle,survival_mean_dle,survival_mean_dfle,
                                      survival_mean_ftc,survival_mean_wo_ftc,survival_pi))
  Outcome<- c("TLE", "DLE", "DFLE","FTC","Wo_FTC","pi")
  survival_mean<- cbind(Outcome, survival_mean)


  writeDataTable(GOM_results, sheet = "Survival", x = survival_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  writeDataTable(GOM_results, sheet = "Survival", x = survival_mean,
                 colNames = TRUE, rowNames = FALSE,startRow = 6+i*6)

  ## transpose survival data to plot average survival curve
  survival_mean_plot<- survival_mean%>%
    gather(time, value, `0`:`10`, factor_key=TRUE)

  output_plot<-ggplot( survival_mean_plot , aes(x=time, y=value, color= Outcome, shape=Outcome,group = Outcome)) +
    geom_point() + geom_line() +
    xlab('Years since Intake') + ylab('Probability of Survival') +
    scale_y_continuous(breaks=seq(0,1,0.1), minor_breaks=NULL, expand=c(0,0)) +
    theme(text = element_text(size=10)) +
    ggtitle("Average survival curve among all people")
  print(output_plot)
  insertPlot(GOM_results, sheet = "Survival",  xy=c("A", 16+i*6))


  #percent
  writeDataTable(GOM_results, sheet = "Percents Table", x = percents_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  #Cond1 Table
  writeDataTable(GOM_results, sheet = "Cond1 Table", x = cond1_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  #Cond2 Table
  writeDataTable(GOM_results, sheet = "Cond2 Table", x = cond2_all,
                 colNames = TRUE, rowNames = FALSE,startRow = 1)

  # save excel
  saveWorkbook(GOM_results, out.excel.file, overwrite = TRUE)

}

