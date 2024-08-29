##RKinetDS v1.2 - software for modeling dissolution profiles
##Authors: 
#Natalia Obajtek: nat.obajtek@gmail.com
#Aleksander Mendyk: mfmendyk@cyf-kr.edu.pl; aleksander.mendyk@uj.edu.pl
#Jakub Szlęk: j.szlek@uj.edu.pl
#Adam Pacławski: adam.paclawski@uj.edu.pl
##License: GPLv3

RKinetDS_comp_core<-function(config_yaml='user_config.yml', progress = NULL){
  
#################################### REQUIRED LIBRARIES ####################################
library(compiler)
library(stringr)
library(yaml)

#################################### CLEAR RESULTS FOLDERS #############################
dirs_to_del <- list.dirs()
dirs_to_del <- grep("Results_eq_*", dirs_to_del, value = TRUE)
unlink(dirs_to_del, recursive=TRUE,force=TRUE)

#################################### LOADING CONFIG ####################################
config <- yaml::yaml.load_file('user_config.yml')

#################################### TEST DATA ####################################
data_file<-read.csv(config$data$input_data, header=config$data$headers,sep="\t", colClasses="numeric")

#################################### OPTIMIZATION  METHODS ####################################
##optim method
use_NM <-config$optim_method$NelderMead
use_SANN <-config$optim_method$SANN
use_nloptr <-config$optim_method$nloptr
use_rgenoud <-config$optim_method$rgenoud
use_gensa <-config$optim_method$genSA

#################################### OPTIMIZATION  OPTIONS ####################################
#BGS_method
maxit_BFGS <-as.numeric(config$running_params$maxit_BFGS)

#Nelder-Mead method
maxit_NM <-as.numeric(config$running_params$maxit_NM)

#SANN_method
maxit_SANN <-as.numeric(config$running_params$maxit_SANN)

#nloptr_method
maxit_nloptr <-as.numeric(config$running_params$maxit_nloptr)
optim_rel_tol <- as.numeric(config$running_params$optim_rel_tol)
opti_trace <- as.numeric(config$running_params$opti_trace)

#rgenoud_method
max_iter_rgenoud <-as.numeric(config$running_params$max_iter_rgenoud)

#gensa_method
max_iter_gensa <-as.numeric(config$running_params$max_iter_gensa)

#################################### EQUATION TO OPTIMIZE ####################################
##Load model equation
selectedModelsIdx <- config$selected_modelsx$selected_models
allModels <- config$models
selectedModels <- allModels[selectedModelsIdx]
selectedModelsNames <- names(selectedModels)
equationToOptim <- unlist(unname(selectedModels))

#################################### LOOP ####################################
max_loop<-1
max_supra_loop<-length(equationToOptim)

#################################### REPORT LISTS ####################################
RMSE_list <- list()
R2_list <- list()
R2_ad_list <- list()
AIC_list <- list()

model_error_list_RMSE <- list()
model_error_list_R2 <- list()
model_error_list_R2_ad <- list()
model_error_list_AIC <- list()

overall_output <- list()

#################################### FUNCTION OF OPTIMIZING ####################################
##RMSE function
RMSE1<-function(matrix, parameters, equat){
  res<-Inf
  C<-parameters
  try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("t", i, sep=""), as.double(matrix[,i]))
    out_RMSE<-as.double(matrix[,dim(matrix)[2]])
  },TRUE)
  
  try(y<-eval(parse(text=equat)),TRUE)
  try(res<-sqrt(mean((y-out_RMSE)^2)),TRUE)
  return(res)
  
}

##Function for minimization
funct1<-function(parameters, equat){
  C<-parameters
  y<-as.numeric(eval(parse(text=equat)))
  res<-mean((y-out)^2)
  if (is.na(res)){res<-Inf}
  return(res)
}

##T-res function
tRes1<-function(matrix, parameters, equat){
  colClasses <- c("numeric", "NULL")
  time <- read.table(config$data$input_data, header=config$data$headers, sep="\t", colClasses=colClasses)
  C<-parameters
  try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("t", i, sep=""), as.double(matrix[,i]))
    observed<-as.double(matrix[,dim(matrix)[2]])
  },TRUE)
  
  try(predicted<-eval(parse(text=equat)),TRUE)
  try(res_t<-cbind(time, observed, predicted), TRUE)
  colnames(res_t)<-c("Time", "Observed", "Predicted")
  return(res_t)
  
}

##R2 function
R_squared <- function(matrix, parameters, equat) {
  res_R2 <- Inf
  C <- parameters
  try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("t", i, sep=""), as.double(matrix[,i]))
    out_R2<-as.double(matrix[,dim(matrix)[2]])
  },TRUE)
  
  y_pred <- eval(parse(text=equat))
  y_obs <- out_R2
  
  ss_residual <- sum((y_obs-y_pred)^2)
  ss_total <- sum((y_obs-mean(y_obs))^2)
  res_R2 <- 1 - (ss_residual/ss_total)
  
  return(res_R2)
  
}

##R2 adjusted function
ad_R_squared <- function(matrix, parameters, equat) {
  res_R2_ad <- Inf
  C <- parameters
  no_observ <- dim(matrix)[1]
  no_param <- length(parameters)
  
  try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("t", i, sep=""), as.double(matrix[,i]))
    out_R2_ad<-matrix[, dim(matrix)[2]]
  },TRUE)
  
  y_pred <- eval(parse(text=equat))
  y_obs <- out_R2_ad
  
  ss_residual <- sum((y_obs-y_pred)^2)
  ss_total <- sum((y_obs-mean(y_obs))^2)
  
  r_squared <- 1 - (ss_residual/ss_total)
  res_R2_ad <- 1-((1-r_squared)*(no_observ-1)/(no_observ-no_param-1))
  
  return(res_R2_ad)
}


##AIC function
AIC1 <- function(matrix, parameters, equat) {
  res_AIC <- Inf
  C <- parameters
  no_observ <- dim(matrix)[1]
  no_param <- length(parameters)
  
  try(for(i in 1:(dim(matrix)[2]-1)) {
    assign(paste("t", i, sep=""), as.double(matrix[,i]))
    out_AIC <- as.double(matrix[,dim(matrix)[2]])
  }, TRUE)
  
  y_pred <- eval(parse(text=equat))
  y_obs <- out_AIC
  ss_residual <- sum((y_obs-y_pred)^2)
  
  res_AIC <- no_observ*log(ss_residual/no_observ)+2*no_param
  return(res_AIC)
}



##Optimize functions
RMSE<-cmpfun(RMSE1)
funct<-cmpfun(funct1)
tRes<-cmpfun(tRes1)
R2<-cmpfun(R_squared)
R2_ad<-cmpfun(ad_R_squared)
AIC<-cmpfun(AIC1)

##Function to assign calculated parameters to equation
get_algebraic_equation<-function(optimizedEquationToSave,C){
  
  for(i in 1:length(C)) {
    pattern<-paste("C[", i, "]", sep="")
    optimizedEquationToSave<-gsub(pattern, C[i], optimizedEquationToSave,fixed = TRUE)
  }
  return(optimizedEquationToSave)
}

############################################################# Preprocesing ################################################################
##Supra optimization loop
for (lk_supra_loop in 1: max_supra_loop) {
  
  single_equationToOptim <-equationToOptim[lk_supra_loop]
  single_equationName <- selectedModelsNames[lk_supra_loop]
  
  ##Correct input equation if you use and 'ln' symbol for logarithm 
  equationBuf<-gsub(" ", "", single_equationToOptim)
  equation<-gsub("ln", "log", equationBuf)
  pattern<-"C\\[\\w+\\]"
  buffer_constants<-unique(unlist(str_extract_all(equation, pattern)))
  N_params<-length(buffer_constants)
  
  ##Prepare parameters and variables
  RMSE_val<-vector(length=max_loop)
  all_params <- matrix(data=NA,nrow=max_loop,ncol=N_params)
  RMSE_ind<-vector(length=max_loop)
  best_all_params<-matrix(data=NA,nrow=max_loop,ncol=N_params)
  best_RMSE_val<-vector(length=max_loop)
  best_RMSE_ind<-vector(length=max_loop)
  best_RMSE_total<-1000000000000
  
  ##Display equation
  print("Equation for optimization")
  print(equation)
  
  ##Main optimization function
  for(lk_loop in 1:max_loop) {
    
    cat("Iteration no = ",lk_loop,"\n")
    cat("Data file =", config$data$input_data, "\n")
    
    assign("t", as.double(data_file[,1]))
    out<-as.double(data_file[,2])
    
    #paramFunct <-vector(length=N_params, "numeric")
    
    paramFunct<-rnorm(N_params)/10
    #paramFunct<-c(8.958148)
    print("paramFunct")    
    print(paramFunct)
    best_error<-100000000
    cat("Iteration no = ",lk_loop,"\n")
    cat("best_error INIT = ",best_error,"\n")
    
    para1.backup<-paramFunct
    X<-data_file
    
    print("Check init values")
    
    preliminary_output<-funct(paramFunct, equation)                      
    cat("Preliminary output = ",preliminary_output,"\n")
    
    #################################### SANNN ####################################
    
    if (use_SANN){
      ## Optim with optim(SANN)
      print(equation)
      print("params on start SANN")
      print(paramFunct)
      
      try(fit0 <- optim(
        paramFunct,
        equat=equation,
        fn=funct,
        method="SANN",
        control=list(trace=opti_trace,maxit=maxit_SANN))
        ,TRUE)
      
      print("summary")
      try(summary(fit0),TRUE)
      
      try(if(fit0$convergence==0){
        par_optim_SANN<-as.vector(as.matrix(fit0$par))
        RMSE_SANN<-RMSE(data_file, par_optim_SANN, equation)
        best_error<-RMSE_SANN
        paramFunct<-par_optim_SANN
        RMSE_val[lk_loop]<-RMSE_SANN
      }else{
        par_optim_SANN<-paramFunct
        RMSE_SANN<-Inf
      },TRUE)
      
    }     
    
    for_domain <- matrix(data=NA,nrow=length(paramFunct),ncol=2)
    
    for(i in 1:length(paramFunct)) {
      for_domain[i,1]<--100*max(abs(paramFunct))
      for_domain[i,2]<-100*max(abs(paramFunct))
    }
    
    #################################### NLOPTR ####################################    
    
    if (use_nloptr){
      
      library(nloptr)
      
      print("Running nloptr")
      fit0 <- nloptr(x0=paramFunct,
                     eval_f=funct,    	
                     lb=for_domain[,1],
                     ub=for_domain[,2],
                     equat=equation,
                     opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
      )
      
      paramFunct<-fit0$solution
      
      print("FINAL RESULTS nloptr")
      print(paramFunct)
      
    }
    
    #################################### GENSA ####################################
    
    if (use_gensa){
      
      library(GenSA)
      
      print("Running GenSA")
      fit0 <- GenSA(par=paramFunct,
                    #errorMeasure=errorOpt,  
                    #ind=inpIndex, 
                    #model=modelOpt, 
                    #tar=tarOpt,    
                    fn=funct,
                    equat=equation,
                    lower=for_domain[,1],
                    upper=for_domain[,2],
                    control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
      )
      
      paramFunct<-fit0$par
      
      print("FINAL RESULTS GenSA")
      print(paramFunct)
      
    }
    
    #################################### RGENOUD ####################################
    
    if (use_rgenoud){
      
      print("Running rgenoud")
      
      library(rgenoud)
      
      fit2 <- genoud(
        fn=funct,	
        nvars=length(paramFunct),
        starting.values=paramFunct,
        equat=equation,
        max.generations=max_iter_rgenoud,
        Domains=for_domain
      )
      
      paramFunct<-fit2$par
      print("Current estimates of the parameters - rgenoud1")
      print(paramFunct)
      
    }
    
    #################################### NELDER-MEAD ####################################
    
    if (use_NM){
      ## Optim with optim(SANN)
      print(equation)
      print("params on start SANN")
      print(paramFunct)
      
      try(fit0 <- optim(
        paramFunct,
        equat=equation,
        fn=funct,
        method="Nelder-Mead",
        control=list(trace=opti_trace,maxit=maxit_NM))
        ,TRUE)
      
      print("summary")
      try(summary(fit0),TRUE)
      
      paramFunct<-fit0$par
      
      print("FINAL RESULTS NM")
      print(paramFunct)
      
    }     
    
    #################################### BFGS ####################################
    
    print("params on start BFGS")
    print(paramFunct)
    
    par_optim_NM<-paramFunct
    
    try(fit1 <- optim(
      paramFunct,
      equat=equation,
      fn=funct,    
      method="BFGS",
      control=list(trace=opti_trace,maxit=maxit_BFGS)
    ),TRUE)
    
    print("FINAL RESULTS OPTIM(BFGS)")
    try(print(fit1$par),TRUE)
    
    
    print("WHOLE object")
    try(print(fit1),TRUE)
    
    try(par_optim_NM<-fit1$par,TRUE)
    
    
    RMSE_NM<-Inf
    print("RMSE")
    try(RMSE_NM<-RMSE(data_file,par_optim_NM, equation))
    print(RMSE_NM)
    
    R2_e<-Inf
    print("R2")
    try(R2_e<-R2(data_file,par_optim_NM, equation))
    print(R2_e)
    
    R2_ad_e<-Inf
    print("R2_adjusted")
    try(R2_ad_e<-R2_ad(data_file,par_optim_NM, equation))
    print(R2_ad_e)
    
    AIC_e<-Inf
    print("AIC")
    try(AIC_e<-AIC(data_file,par_optim_NM, equation))
    print(AIC_e)
    
    cat("Iteration no = ",lk_loop,"\n")
    print("Final params")
    try(print(par_optim_NM),TRUE)
    try(all_params[lk_loop,]<-par_optim_NM,TRUE)
    print(" ")
    
    ##This conditions should help with situation where test RMSE is NA what cause problems and stop current job.
    if(is.na(RMSE_NM)){
      RMSE_NM<-Inf
    }
    
    if(is.na(R2_e)){
      R2_e<-Inf
    }
    
    if(is.na(R2_ad_e)){
      R2_ad_e<-Inf
    }
    
    if(is.na(AIC_e)){
      AIC_e<-Inf
    }
    
    rmserror<-RMSE_NM
    RMSE_ind[lk_loop]<-rmserror
    
    R2error<-R2_e
    RMSE_ind[lk_loop]<-R2error
    
    R2aderror<-R2_ad_e
    RMSE_ind[lk_loop]<-R2aderror
    
    AICerror<-AIC_e
    RMSE_ind[lk_loop]<-AICerror
 
    cat("Iteration no = ",lk_loop,"\n")
    cat("RMSE_test = ",rmserror,"\n")
    cat("R2_test = ",R2error,"\n")
    cat("R2_adjusted_test = ",R2aderror,"\n")
    cat("AIC_test = ",AICerror,"\n")
    
    # try(save(fit1,file=outfile),TRUE)
    
    print("-------------------------------------")
    
  }         
  ##End of optimization function
  ##----------------------------------------------------
  
  ## End of lk_loop

  ##Assign calculated parameters to equation
  algebraic_equation<-get_algebraic_equation(single_equationToOptim,paramFunct) 
  
  #################################### SUMMARY ####################################
  
  cat("\n\n", "SUMMARY", "\n")
  cat("Model name: ", single_equationName, "\n")
  cat("Equation: ", single_equationToOptim, "\n")
  cat("Params:", paramFunct, sep="\t", "\n")
  cat("Algebraic form: ", algebraic_equation, "\n")
  cat("RMSE: ", rmserror, "\n")
  cat("R2: ", R2error, "\n")
  cat("R2 adjusted: ", R2aderror, "\n")
  cat("AIC: ", AICerror, "\n")
  
  ##Save equation and params into text files
  
  folder_names<-c(paste("Results_eq_", selectedModelsIdx[lk_supra_loop],sep=""))
  folder_result<-(file.path(folder_names))
  if (dir.exists(folder_result)) {
    unlink(folder_result, recursive=TRUE,force=TRUE)
    dir.create(folder_result)
  } else {
    dir.create(folder_result)
  }
  
  name_directories<-paste(folder_result, "/optimizedEquation_", selectedModelsIdx[lk_supra_loop], ".txt", sep="")
  sink(as.character(name_directories))
  cat("Model name: ", single_equationName, "\n\n")
  cat("Model equation: ", single_equationToOptim, "\n\n")
  cat("Parameters for equation after optimization: ")
  cat(paramFunct, sep="\t")
  cat("\n\n")
  cat("Algebraic form: ", algebraic_equation, "\n\n")
  cat("RMSE: ", rmserror, "\n\n")
  cat("R2:", R2error, "\n\n")
  cat("R2 adjusted:", R2aderror, "\n\n")
  cat("AIC:", AICerror)
  cat("\n\n")
  sink()
  
  obsPredFileName<-paste(folder_result, "/results_", selectedModelsIdx[lk_supra_loop], ".txt", sep="")
  predObsMat<-tRes(
    matrix=data_file,
    parameters=paramFunct,
    equat=equation)
  col.names <-c("Observed", "Predicted")
  
  sink(as.character(obsPredFileName))
  write.table(predObsMat, sep="\t", file=obsPredFileName, col.names=TRUE, row.names=FALSE, append=FALSE)
  sink()

  ##Error and report list
  #RMSE
  model_error_single_RMSE <- paste(selectedModelsIdx[lk_supra_loop],single_equationToOptim,rmserror,sep="\t")
  model_error_list_RMSE[lk_supra_loop] <- paste(model_error_single_RMSE,"\n")
  RMSE_list[lk_supra_loop] <- c(rmserror)
  
  #R2
  model_error_single_R2 <- paste(selectedModelsIdx[lk_supra_loop],single_equationToOptim,R2error,sep="\t")
  model_error_list_R2[lk_supra_loop] <- paste(model_error_single_R2,"\n")
  R2_list[lk_supra_loop] <- c(R2error)
  
  #R2 adjusted
  model_error_single_R2_ad <- paste(selectedModelsIdx[lk_supra_loop],single_equationToOptim,R2aderror,sep="\t")
  model_error_list_R2_ad[lk_supra_loop] <- paste(model_error_single_R2_ad,"\n")
  R2_ad_list[lk_supra_loop] <- c(R2aderror)
  
  #AIC
  model_error_single_AIC <- paste(selectedModelsIdx[lk_supra_loop],single_equationToOptim,AICerror,sep="\t")
  model_error_list_AIC[lk_supra_loop] <- paste(model_error_single_AIC,"\n")
  AIC_list[lk_supra_loop] <- c(AICerror)  
  
  overall_output[[lk_supra_loop]] <- list(
  "Model name: ",single_equationName,"\n","Model equation:",single_equationToOptim,"\n",
  "Parameters for equation after optimization: ",paramFunct,"\n",
  "Algebraic form: ", algebraic_equation, "\n",
  "RMSE: ", rmserror,"\n",
  "R2: ", R2error,"\n",
  "R2 adjusted: ", R2aderror,"\n",
  "AIC: ", AICerror, "\n\n")
  
  # Updating progress bar - shiny only function
  # If we passed a progress update function, call it
  try(if (progress) {
    incProgress(amount = 1/length(equationToOptim),
                detail = paste0('Calculating eq: ', lk_supra_loop, ' out of ',
                  length(equationToOptim)))
  }, TRUE)

}
## End of loop lk_supra_loop

##Create a ranking based on RMSE
RMSE_order <- order(unlist(RMSE_list))
sorted_data_RMSE <- model_error_list_RMSE[RMSE_order]

##Create a ranking based on R2
R2_order <- order(unlist(R2_list), decreasing = TRUE)
sorted_data_R2 <- model_error_list_R2[R2_order]

##Create a ranking based on R2 adjusted
R2_ad_order <- order(unlist(R2_ad_list), decreasing = TRUE)
sorted_data_R2_ad <- model_error_list_R2_ad[R2_ad_order]

##Create a ranking based on AIC
AIC_order <- order(unlist(AIC_list))
sorted_data_AIC <- model_error_list_AIC[AIC_order]

#Report RMSE ranking
sink("RMSE_error_ranking.txt")
cat("Number\tEquation\tRMSE\n")
cat(unlist(sorted_data_RMSE),sep="")
sink()

#Report R2 ranking
sink("R2_error_ranking.txt")
cat("Number\tEquation\tR2\n")
cat(unlist(sorted_data_R2),sep="")
sink()

#Report R2 adjusted ranking
sink("R2_adjusted_error_ranking.txt")
cat("Number\tEquation\tR2_adjusted\n")
cat(unlist(sorted_data_R2_ad),sep="")
sink()

#Report AIC ranking
sink("AIC_error_ranking.txt")
cat("Number\tEquation\tAIC\n")
cat(unlist(sorted_data_AIC),sep="")
sink()

##Make overall report
sink("overall_report.txt")
cat("OVERALL REPORT", "\n\n")
cat("Date: ", date(), "\n\n")
cat("USED OPTIMIZATION METHODS", "\n")
cat("NelderMead: ", config$optim_method$NelderMead, "\n")
cat("SANN: ", config$optim_method$SANN, "\n")
cat("nloptr: ", config$optim_method$nloptr, "\n")
cat("rgenoud: ", config$optim_method$rgenoud, "\n")
cat("genSA: ", config$optim_method$genSA, "\n\n")
cat("USED PARAMETERES", "\n")
cat("maxit_BFGS: ", maxit_BFGS, "\n")
cat("maxit_NM: ", maxit_NM, "\n")
cat("maxit_nloptr: ", maxit_nloptr, "\n")
cat("max_iter_rgenoud: ", max_iter_rgenoud, "\n")
cat("max_iter_gensa: ", max_iter_gensa, "\n")
cat("maxit_SANN: ", maxit_SANN, "\n")
cat("optim_rel_tol: ", optim_rel_tol, "\n")
cat("opti_trace: ", opti_trace, "\n\n")
cat("RESULTS", "\n")
cat(unlist(overall_output),sep="")
cat("","\n")
cat("SYSTEM USER INFO", "\n")
cat(paste(Sys.info()), sep="\t")
cat("","\n")
cat("","\n")
cat("INSTALLED PACKAGES AND THEIR VERSIONS ","\n")
ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]

for(i_ip in 1:nrow(ip)){
  cat(unlist(ip[i_ip,]),sep="\t","\n")
}
sink()
}
