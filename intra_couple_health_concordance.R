#INFO
########################################################################################################
#
# Intra-couple health concordances among Hungarian couples
# R-Code: Zoltan Brys, Gergely Tóth, Gergely Rosta
#
#
########################################################################################################


# ENVIRONMENT, PACKAGES, FUNCTIONS
########################################################################################################
#environment
  rm(list = ls())
  if (as.numeric(gsub(".*:(\\s*)(\\d+)(\\s+)\\d+.*", "\\2", (system("free -m", intern = TRUE)[2])))<1024) 
    stop("Memory are likely not will be enough!")
  if (!("stats" %in% (.packages()) )) stop("R Environment is not fully loaded!")

#packages
  library("gtsummary")   #descriptive tables

#raw data URL from KDK data repository
  raw_data_url <- url("https://openarchive.tk.mta.hu/595/5/2023_dataset.RDS", "rb")

#data manipulation function recode factors into binary with NAS
recode_bin <- function(f, ones, zeros, nas) 
  {
    result <- rep(NA, length(f))  # Initialize result vector with NA
    
    for (i in 1:length(f)) {
      if (f[i] %in% zeros) {
        result[i] <- 0
      } else if (f[i] %in% ones) {
        result[i] <- 1
      } else if (f[i] %in% nas) {
        result[i] <- NA
      }
    }
      return(result)
} #end of recode_bin


#one sample proportion 99% CI estimation with z-test
 #p is the probability
 #n is the sample size
one_z <- function (p, n, z=2.576)  
{
  
  if (length(p)!=1) stop("p is vector!")
  if (length(n)!=1) stop("n is vector!")
  
  p1    = p                  # sample p
  pq    = p1 * (1-p1)        # sample p*(1-p)
  pqn   = pq / n             # sample p*(1-p) / n
  spqn  = sqrt(pqn)          # sample sqrt () p*(1-p) / n)
  zspqn = spqn *  z          #99% CI with default z
  
  df_res = data.frame(p = p1 , 
                      pl = (p1 - zspqn) ,
                      ph = (p1 + zspqn))
  
  return(df_res) #p, low , high
} #end of one sample z-test


#one sample z test based on binary vector
 #x is binary vector
 #nas = if NAs allowed
 #z= CI% z-value, default z is 99% CI
one_z_binvec <- function (x, nas = TRUE, z=2.576)  
{
  
  if (nas == FALSE)
  {
    if (anyNA(x)) stop("NAs present in x!")
  }
  
  if (nas == TRUE) 
  { 
    dfx <- data.frame(x=x)
    dfx <- na.omit(dfx)
    x <- dfx$x
  }
  
  p1 = sum(x) / length(x) # sample p
  df_res = one_z(p1, length(x), z=z) 
  
  return (df_res) #p, low , high
} #end of one sample z-test

#two sample z test and Cohen's h based on p values
  #px is prob of first sample, nx = size of first sample
  #py is prob of second sample, ny = size of second sample
  #z= CI% z-value, default z is 99% CI
two_zh <- function (px, py, nx, ny, z=2.576)   
{
  
  if (length(px)!=1) stop("px is vector!")
  if (length(py)!=1) stop("py is vector!")
  if (length(nx)!=1) stop("nx is vector!")
  if (length(ny)!=1) stop("ny is vector!")
  
  dp = px - py #dif of tow sample ps
  zsedp = (sqrt(px*(1-px)/nx+py*(1-py)/ny))*z
  
  #Cohen's H
  ch = ( (2*asin(sqrt(px))) - (2*asin(sqrt(py))) )
  se = sqrt(0.25 * (1 / nx + 1 / ny ))
  mo = z * se

  df_res = data.frame(dp=dp, 
                      dpl = (dp - zsedp), 
                      dph = (dp + zsedp), 
                      ch  = ch,
                      chl =  (ch - mo),
                      chh =  (ch + mo)
                      )
  
  return (df_res) #dp, dp low, dp high, Cohen's h for proportion
  #Cohen's h for proportion ranges -pi +pi!
} #end of two sample z test and Cohen's h for proportion


#two sample z test and Cohen's h based on binary vectors
  #x binary vector
  #y binary vector
  #nas = if NAs allowed (omitted)
  #z= CI% z-value, default z is 99% CI
two_zh_bin_vecs <- function (x, y, nas = TRUE, z=2.576)   
{
  
  if (length(x) != length(y)) stop("Non-equal lenghts.")
  
  if (nas == FALSE)
  {
    if (anyNA(x)) stop("NAs present in x!")
    if (anyNA(y)) stop("NAs present in y!")
  }
  
  
  if (nas == TRUE) 
  { 
    dfxy <- data.frame(x=x, y=y)
    dfxy <- na.omit(dfxy)
    x <- dfxy$x
    y <- dfxy$y
  }
  
  nx = length(x) #sample size
  ny = length(y) #sample size , this case it is equal
  
  px = sum(x) / nx # sample p x
  py = sum(y) / ny # sample p y
  
  df_res = two_zh(px, py, nx, ny, z=z)
  
  return(df_res)
 #Cohen's h for proportion ranges -pi +pi!
} #end of two sample z test and Cohen's h for proportion


#main function concordance measures
 #independent probability for concordance - random matching and 99% CI
 #obsereved probability of concordance and99% CI
 #difference between the previous two ps with 99% CI
 #Cohen's h between the prev two (for proportion ranges between -pi and +pi)
conc_mes <- function (x, y, nas = TRUE, z=2.576)   
{
  
  if (length(x) != length(y)) stop("Non-equal lenghts.")
  
  if (nas == FALSE)
  {
    if (anyNA(x)) stop("NAs present in x!")
    if (anyNA(y)) stop("NAs present in y!")
  }
  
  if (nas == TRUE) 
  { 
    dfxy <- data.frame(x=x, y=y)
    dfxy <- na.omit(dfxy)
    x <- dfxy$x
    y <- dfxy$y
  }
  
  n= length(x)

  #independent probability for concordance - random matching
  pivc = (sum(x==1) / n) * (sum(y==1) / n) + (sum(x==0) / n) * (sum(y==0) / n)
  
  tmp_ip_c = one_z(pivc, (2*n) ) #99% CI for probability of concordance w def z
  
  icp  =  tmp_ip_c$p  #indepedent concordance p 
  icpl =  tmp_ip_c$pl #indepedent concordance p low
  icph =  tmp_ip_c$ph #indepedent concordance p high

  #obsereved concordance
  c11 = sum(x == 1 & y == 1)
  c00 = sum(x == 0 & y == 0)

  op = (c11+c00) / n
  tmp_op_c = one_z(op, n ) #99% CI of observed concordance
  
  ocp  = tmp_op_c$p  #observerd concordance p 
  ocpl = tmp_op_c$pl #observerd concordance p low
  ocph = tmp_op_c$ph #observerd concordance p high
  
  #difference 
  tmp_two_zh = two_zh(ocp, icp, n, n)
  dp  =  tmp_two_zh$dp
  dpl =  tmp_two_zh$dpl
  dph =  tmp_two_zh$dph 
  
  #cohen's h
  ch  =  tmp_two_zh$ch
  chl =  tmp_two_zh$chl
  chh =  tmp_two_zh$chh
  
  
  res_df = data.frame(icp, icpl, icph, ocp, ocpl, ocph, dp, dpl, dph, ch, chl, chh)
    #icp = independent concordance
    #ocp = observed concordance
    #dp = difference
    #Cohen's
  return(res_df)
}
########################################################################################################


# READ RAW DATA FROM PUBLIC SOURCE
########################################################################################################
#read RDS
  survey_raw <- readRDS(raw_data_url)
  if (!exists("survey_raw")) stop("Input data RDS from KDK could not be read!")

#convert to data frame
  df_survey <- as.data.frame(survey_raw)
  if (!exists("df_survey")) stop("Input data is incorrect!")

#check the rows/columns
  if (dim(df_survey)[1] != 1500) stop("Input data is incorrect!")
  if (dim(df_survey)[2] != 195) stop("Input data is incorrect!")

#close connection and delete temporary vars
  close(raw_data_url)
  rm(survey_raw, raw_data_url)
########################################################################################################

  
# SELECT VARIABLES OF INTEREST
########################################################################################################
#variables of interest
  at_var_names <-  c("szulev", #birth year
                     "SZD1",   #gender
                     "SZD2",   #number of adults
                     "SZD2g",  #number of minors
                     "SZD4",   #education
                     "SZD9",   #financial status
                     "SZD6a",  #settlement type
                     "SZD6b",  #county
                     
                     "par",    #couple status, 1=living together
                     "P1_ev",  #length of relationship (in years)
                     "P2",     #partner : birth of year
                     "P3",     #partner educational level
                     
                     "HP3",    #partner perceived smoking
                     "HP4",    #partner perceived influenza-occulation
                     "HP7",    #partner perceived vegetable
                     "HP8",    #partner perceived sleeping
                     "HP9",    #partner perceived exercise
                     "pMV1",   #partner perceived covid19
                     "pSZ1",   #partner perceived future hypotethical vacc
                     
                     
                     "H3",    #respondent smoking
                     "H4",    #respondent influenza-occulation
                     "H7",    #respondent vegetable
                     "H8",    #respondent sleep
                     "H9",   #respondent exercise
                     "MV1",    #respondent covid19
                     "SZ1"    #respondent future hypotethical vacc
                    )   
                  
#subsetting to variables of interest 
  df <- subset(df_survey, select = at_var_names)

#delete unused vars
  rm(df_survey)
########################################################################################################  


# SELECT PARTICIPANTS
#######################################################################################################  
#select participants
  df$filter1 <- (df$par == "Igen, együtt élünk.") & ( df$SZD2==2 )
  table(df$filter1) #prints number of participants under TRUE
  pc <- subset(df, subset = (filter1))
  rm(df) #delete full
####################################################################################################### 
 
  
# DATA PREPARATION 
#######################################################################################################  
#handle NAs of factors
cols_fac <- colnames(pc)[sapply(pc, is.factor)] #factor colnames

for (col in cols_fac) {
  summary(pc[[col]])
  pc[[col]][pc[[col]] == "Nem tudom" | pc[[col]] == "Nem akarok válaszolni"] <- NA
  pc[[col]] <- droplevels(pc[[col]], exclude = c("Nem tudom", "Nem akarok válaszolni"))
}

#binarized tobacco use respondent and partner
  #levels(pc$H3)
  pc$rhb1 <- recode_bin( f =  pc$H3,
              ones = c("Dohányzom, de elektronikus eszközt (e-cigaretta, hevített dohánytermék) nem használok.", 
                       "Dohányzom és elektronikus eszközt (e-cigaretta, hevített dohánytermék) is használok.", 
                       "Csak elektronikus eszközt (e-cigaretta, hevített dohánytermék) használok.") , 
              zeros = c("Már leszoktam a dohányzásról, és/vagy az elektronikus eszköz (e-cigaretta, hevített dohánytermék).",
                        "Sohasem dohányoztam és elektronikus eszközt (e-cigaretta, hevített dohánytermék) sem használtam."),
              nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
             )
    attr(pc$rhb1, "label") <- "Tobacco use respondent (0-1)"

  #levels(pc$HP3)
  pc$phb1 <- recode_bin( f =  pc$HP3,
                         ones = c("Dohányzik, de elektronikus eszközt nem használ.", 
                                  "Dohányzik és elektronikus eszközt is használ.", 
                                  "Csak elektronikus eszközt használ.") , 
                         zeros = c("Már leszokott a dohányzásról, és/vagy az elektronikus eszköz használatról.",
                                   "Sohasem dohányozott és elektronikus eszközt sem használt soha."),
                         nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                        )
    attr(pc$phb1, "label") <- "Tobacco use partner (0-1)"

#binarized covid-19 binarized respondent and partner
  #levels(pc$MV1)
  pc$rhb2 <- recode_bin( f =  pc$MV1,
                        ones = c("Nem és nem is tervezem. ", 
                                 "Még nem, de tervezem. ") , 
                        zeros = c("Igen, és tervezek továbbiakat is. ",
                                  "Igen, de nem tervezek továbbiakat. "),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                        )
    attr(pc$rhb2, "label") <- "Covid-19 vaccination hesitancy respondent (0-1)"

  #levels(pc$pMV1)
  pc$phb2 <- recode_bin( f =  pc$pMV1,
                         ones = c("Nem és nem is tervezi.", 
                                  "Még nem, de tervezi.") , 
                         zeros = c("Igen, és tervez továbbiakat is.",
                                   "Igen, de nem tervez továbbiakat."),
                         nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                        )
    attr(pc$phb2, "label") <- "Covid-19 vaccination hesitancy partner (0-1)"

#flu-inocculation
  #levels(pc$H4)  
  pc$rhb3 <- recode_bin( f =  pc$H4,
                         ones = c("Nem, még soha nem adattam be magamnak.") , 
                         zeros = c("Igen, rendszeresen, minden évben.",
                                   "Igen, időnként.",
                                   "Egyszer vagy kétszer beadattam."),
                         nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                        )
    attr(pc$rhb3, "label") <- "Influenza vaccination uptake never respondent (0-1)"

  #levels(pc$HP4)  
  pc$phb3 <- recode_bin( f =  pc$HP4,
                         ones = c("Nem, még soha nem adatta be magának.") , 
                         zeros = c("Igen, rendszeresen, minden évben.",
                                   "Igen, időnként.",
                                   "Egyszer vagy kétszer beadatta."),
                         nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                        )
    attr(pc$phb3, "label") <- "Inluenza vaccination uptake never partner (0-1)"
  
#vegetable
 #levels(pc$H7)
 pc$rhb4 <- recode_bin( f =  pc$H7,
                        ones = c("hetente",
                                 "ritkábban") , 
                        zeros = c("naponta többször",
                                  "naponta",
                                  "hetente többször"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                         )
  attr(pc$rhb4, "label") <- "Not enough vegetable consumption respondent (0-1)"
 
 #levels(pc$HP7)
 pc$phb4 <- recode_bin( f =  pc$HP7,
                        ones = c("hetente",
                                 "ritkábban") , 
                        zeros = c("naponta többször",
                                  "naponta",
                                  "hetente többször"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
                       )
  attr(pc$phb4, "label") <- "Not enough vegetable consumption partner (0-1)"

#physical exercise
 #levels(pc$H9)
 pc$rhb5 <- recode_bin( f =  pc$H9,
                        ones = c("soha",
                                 "ritkábban, mint heti egyszer"),
                        zeros = c("hetente egyszer",
                                  "hetente többször",
                                  "naponta egyszer",
                                  "naponta többször",
                                  "naponta"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
  attr(pc$rhb5, "label") <- "Not enough physical excercise respondent (0-1)"
 
 #levels(pc$HP9)
 pc$phb5 <- recode_bin( f =  pc$HP9,
                        ones = c("soha",
                                 "ritkábban, mint heti egyszer"),
                        zeros = c("hetente egyszer",
                                  "hetente többször",
                                  "naponta egyszer",
                                  "naponta többször",
                                  "naponta"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
   attr(pc$phb5, "label") <- "Not enough physical excercise partner (0-1)"

#sleep
 #levels(pc$H8)
 pc$rhb6 <- recode_bin( f =  pc$H8,
                        ones = c("Hetente legalább háromszor",
                                 "Hetente"),
                        zeros = c("Havonta vagy ritkábban",
                                  "Soha"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
   attr(pc$rhb6, "label") <- "Regular sleep problems respondent (1-4)"
 
 #levels(pc$HP8)
 pc$phb6 <- recode_bin( f =  pc$HP8,
                        ones = c("Hetente legalább háromszor",
                                 "Hetente"),
                        zeros = c("Havonta vagy ritkábban",
                                  "Soha"),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
    attr(pc$phb6, "label") <- "Regular sleep problems partner (1-4)"

#hypothethical vaccination
 #levels(pc$SZ1)
 pc$rhb7 <- recode_bin( f =  pc$SZ1,
                        ones = c("Biztosan nem oltatnám be magam. ",
                                 "Valószínűleg nem oltatnám be magam. "),
                        zeros = c("Biztosan beoltatnám magam. ",
                                  "Valószínűleg beoltatnám magam. "),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
   attr(pc$rhb7, "label") <- "Hypothetical vaccination hesitency respondent (0-1)"

 #levels(pc$pSZ1)
 pc$phb7 <- recode_bin( f =  pc$pSZ1,
                        ones = c("Biztosan nem oltatná be magát.",
                                 "Valószínűleg nem oltatná be magát."),
                        zeros = c("Biztosan beoltatná magát." ,
                                  "Valószínűleg beoltatná magát."),
                        nas = c(NA, NaN, "NA", "Nem tudom.", "Nem akarok válaszolni.")
 )
    attr(pc$phb7, "label") <- "Hypothetical vaccination hesitency partner (0-1)"
####################################################################################################### 

 
# DESCRIPTIVE VARIABLES
#######################################################################################################
# x1 = sum age of the couple
  pc$x1 <- (2023 - pc$szulev) + (2023- pc$P2)
  attr(pc$x1, "label") <- "Sum ages of the couple (in years)"
# x2 = difference of ages
  pc$x2 <- abs((2023 - pc$szulev) - (2023- pc$P2))
  attr(pc$x2, "label") <- "Difference of partners' ages (in years)"
# x3 = sum edu level
  pc$P3 <- as.numeric(pc$P3)
  pc$P3[pc$P3==8] <- NA #NAs
  pc$x3 <- as.numeric(pc$SZD4) + pc$P3
  attr(pc$x3, "label") <- "Sum education of the couple (2-14)"
# x3 = edu level difference
  pc$x4 <- abs(as.numeric(pc$SZD4) - pc$P3)
  attr(pc$x4, "label") <- "Difference of education of the partners (0-6)"
# number of minors living with the couple
  pc$x5 <- pc$SZD2g
  attr(pc$x5, "label") <- "Number of minors living with the couple (0- )"
# financial status
  pc$x6 <- as.numeric(pc$SZD9)
  pc$x6[pc$x6==6]<- NA
  pc$x6 <- 6 - pc$x6 
  attr(pc$x6, "label") <- "Financial status of the couple (1-5)" 
    #5 = "gondok nélkül él/élnek, "
    #4 = "beosztással jól kijön/kijönnek, "                    
    #3 = "éppen hogy kijön/kijönnek a havi jövedelmükből, "
    #2 = "hónapról-hónapra anyagi gondjai(k) van/vannak, vagy "
    #1 = "nélkülözések között él/élnek? "
    
# settlement type
  pc$x7 <- 5 - as.numeric(pc$SZD6a)
  attr(pc$x7, "label") <- "Settlement type/size (1-4)"
    # 4 = "Budapest "      
    # 3 = "megyeszékhely " 
    # 2 = "város "
    # 1 = "falu " 
#time of cohab
  pc$x8 <- pc$P1_ev
  attr(pc$x8, "label") <- "Cohabiting duration (in years)"

#meta
  xvars <- c("x1","x2","x3","x4","x5","x6","x7","x8")
  xlabels <- c("Sum ages of the couple (in years)",
               "Difference of partners' ages (in years)",
               "Sum education of the couple (2-14)",
               "Difference of education of the partners (0-6)",
               "Number of minors living with the couple (0- )",
               "Financial status of the couple (1-5)",
               "Settlement type/size (1-4)",
               "Cohabiting duration (in years)"
               )
####################################################################################################### 

# TABLE1 - Characteristics of the sample
#######################################################################################################
  
#labels defined for tbl_sum
  flabels <- list(x1 ~ "Sum ages of the couple (in years)",
                  x2 ~ "Difference of partners' ages (in years)",
                  x3 ~ "Sum education of the couple (2-14)",
                  x4 ~ "Difference of education of the partners (0-6)",
                  x5 ~ "Number of minors living with the couple (0- )",
                  x6 ~ "Financial status of the couple (1-5)",
                  x7 ~ "Settlement type/size (1-4)",
                  x8 ~ "Cohabiting duration (in years)"
                  )
#how to report
  fstat <- list(all_continuous() ~ "{median} [{p25} — {p75}]")
    
  tbl_summary(pc,
              label = flabels,
              include = all_of(xvars),
              statistic = fstat,
              digits = list(all_continuous() ~ c(0, 0))
              )
#######################################################################################################


# TABLE2 - concordance if independent, observed concordance, z-test, Cohen's H
#######################################################################################################

#define table2
  table2 <- rbind( data.frame(var="hb1",label="Tobacco use",conc_mes(pc$rhb1, pc$phb1) ) ,
                  data.frame(var="hb2",label="Covid-19 vaccination uptake",conc_mes(pc$rhb2, pc$phb2) ) ,
                  data.frame(var="hb3",label="Flu vaccination uptake",conc_mes(pc$rhb3, pc$phb3) ) , 
                  data.frame(var="hb4",label="Vegetable consumption",conc_mes(pc$rhb4, pc$phb4) ) , 
                  data.frame(var="hb5",label="Physical excercise",conc_mes(pc$rhb5, pc$phb5) ) , 
                  data.frame(var="hb6",label="Sleep problems",conc_mes(pc$rhb6, pc$phb6) ) , 
                  data.frame(var="hb7",label="Hypothetical vaccination \n uptake intetnion",conc_mes(pc$rhb7, pc$phb7) )
                 )

#sort by desceding ch
  table2 <- table2[order(table2$ch, decreasing = TRUE),]
  
#adding attributes   
  attr(table2$icp, "label")   <- "Concordance if matching is independent"
  attr(table2$icpl, "label")  <- "Concordance if matching is independent 99% CI low"
  attr(table2$icph, "label")  <- "Concordance if matching is independent 99% CI high"
  attr(table2$ocp, "label")   <- "Observed concordance"
  attr(table2$ocpl, "label")  <- "Observed concordance 99% CI low"
  attr(table2$ocph, "label")  <- "Observed concordance 99% CI high"
  attr(table2$dp, "label")    <- "Difference"
  attr(table2$dpl, "label")   <- "Difference 99% CI low"
  attr(table2$dph, "label")   <- "Difference 99% CI high"
  attr(table2$ch, "label")    <- "Cohen's H "
  attr(table2$chl, "label")    <- "Cohen's H 99% low"
  attr(table2$chh, "label")    <- "Cohen's H 99% high"
  
  print(table2, digits = 2)
#######################################################################################################


#Figure 1.    
#######################################################################################################  
library(ggplot2)
  
plot1 <- ggplot(table2, aes(x = reorder(label, -ch)))  +
  geom_point(aes(y = ch), color = "blue",  size = 2) +
  geom_point(aes(y = chl), color = "red") +
  geom_point(aes(y = chh), color = "red") +
  geom_segment(aes(x = reorder(label, -ch), xend = reorder(label, -ch), y = chl, yend = chh), color = "black") +
  labs(x = "Health domain", y = "Cohen's H with 99% CI") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,  size = 12))
  
ggsave("plot1.png", plot = plot1)
#######################################################################################################
  
