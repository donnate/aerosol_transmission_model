##### Simplified COVID-risk estimator
RISKS <-  read_csv("CDC_risk_age.csv",col_types = cols(
  `Age Group` = col_character(),
  Risk = col_double(),
  Type = col_character(),
  Effect = col_character()
))


AGE_CAT= unique(RISKS$`Age Group`)
hospitalization_risk <- function(age){
  cat = 1 + (age > 4) + (age > 17) + (age > 29) + (age > 39)  +
        (age > 49) + (age > 64) + (age > 74) + (age  > 84)
  a = RISKS %>% filter(`Age Group` == AGE_CAT[cat], Type == "hospitalization")
  b = RISKS %>% filter(`Age Group` == AGE_CAT[cat], Type == "death")
  risk_hosp = 17.2
  risk_death = 0.2
  
  risk_hosp = ifelse(a$Effect == "higher",  (risk_hosp * a$Risk),
     ifelse(a$Effect == "lower", risk_hosp / a$Risk, risk_hosp))
  risk_death = ifelse( (b$Effect == "higher"), (risk_death * b$Risk), 
    ifelse((b$Effect == "lower") , (risk_death / b$Risk),risk_death))
   return(c(risk_hosp, risk_death)/100)
  
}