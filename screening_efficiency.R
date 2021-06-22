###### Computes the probability of failing the test

source("proba_symptoms.R")
source("beta_params.R")
### p_asypt = https://pubmed.ncbi.nlm.nih.gov/32691881/ 15.6% (95% CI, 10.1%-23.0%)

data_efficiency_LFA_antigen_asymptomatic <- data.frame("average"= c(80, 40),
                                                       "worst_case" = c(33,18), 
                                                       "best_case" = c(91, 67), 
                                                       row.names = c("Sympt", "Asympt")) *0.01
data_efficiency_LFA_antigen  <- read_csv("lfa_sensitivity.txt")
PROBA_SYMPTOMS = compute_proba_symptoms(lag=input$time2event)
beta_p_asympt_params = beta.parms.from.quantiles(c(0.101, 
                                                 0.23))
beta_p_asympt_test_params = beta.parms.from.quantiles(c(0.18, 
                                                   0.67))
 
compute_sensitivity <- function(P_LIE, input, period4predicting, B=200){
  PERIOD_FOR_PREDICTING = period4predicting
  SENSITIVITY = sapply(1:B, function(b){
    c(sapply(1: (PERIOD_FOR_PREDICTING - input$time2event -1), function(i){
      #print(i)
      time2test = min(PERIOD_FOR_PREDICTING -input$time2event - i, length(PROBA_SYMPTOMS))
      if (time2test>2){
        beta_sympt_params = beta.parms.from.quantiles(c(1-0.01*data_efficiency_LFA_antigen$upper[min(21, time2test)], 
                                                        1-0.01*data_efficiency_LFA_antigen$lower[min(21,time2test)]))
        p_test_sympt =  rbeta(1, beta_sympt_params$a, beta_sympt_params$b)
      }else{
        p_test_sympt =  runif(1, 1-0.01*data_efficiency_LFA_antigen$upper[min(21,time2test)], 1-0.01*data_efficiency_LFA_antigen$lower[min(21,time2test)])
      }
      
      p_asympt =   rbeta(1, beta_p_asympt_params$a, beta_p_asympt_params$b)
      
      p_test_asympt =  rbeta(1, beta_p_asympt_test_params$a, beta_p_asympt_test_params$b)
      av = (1-p_asympt) *   (PROBA_SYMPTOMS[time2test + 2]* P_LIE + (1-PROBA_SYMPTOMS[time2test + 2])) *  (1-p_test_sympt)  +
            p_asympt  * (1-p_test_asympt) 
      return(av)
    }),
    sapply((PERIOD_FOR_PREDICTING - input$time2event): (PERIOD_FOR_PREDICTING-1), function(i){
      time2event = PERIOD_FOR_PREDICTING - i 
      p_asympt =   rbeta(1, beta_p_asympt_params$a, beta_p_asympt_params$b)
      av = (1-p_asympt) *   PROBA_SYMPTOMS[time2event] *  P_LIE +
           ((1-p_asympt) *   (1- PROBA_SYMPTOMS[time2event]) + p_asympt)
      return(av)
    }), 
    1.0
    )
    })
  
  SENSITIVITY_MEAN = apply( SENSITIVITY, 1, mean)
  SENSITIVITY_SD = apply( SENSITIVITY, 1, sd)
  SENSITIVITY_q025 = apply( SENSITIVITY, 1, quantile, 0.025, na.rm=TRUE)
  SENSITIVITY_q975 = apply( SENSITIVITY, 1, quantile, 0.975, na.rm=TRUE)
  return(data.frame("sens"=SENSITIVITY_MEAN,  "q025"= SENSITIVITY_q025,
              "q975"= SENSITIVITY_q975,
              "sd"= SENSITIVITY_SD))
}


