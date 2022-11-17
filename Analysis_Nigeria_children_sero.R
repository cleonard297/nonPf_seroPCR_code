#******************************************************************#
#      ANALYSIS OF NIGERIA CHILDRENS MALARIA SEROLOGICAL DATA      #
#******************************************************************#
# author: Colleen Leonard
# date: November 11, 2021

#Load required packages
library(haven)
library(survey) #to analyze complex survey data
library(srvyr) #dplyr-like syntax for summary stats of survey data
library(plyr) 
library(rgdal)
library(lme4)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(effects)
library(ggeffects)
library(gridExtra)
library(readxl)
library(janitor)
library(skimr)
library(expss)
library(performance)
library(writexl)
library(testthat)
library(aods3) #checks for goodness of fit of GLMM
library(parameters) #checks model parameters
library(scales)
library(labelled)
library(readxl)
library(knitr)
                        
# get working drive
getwd() 

# set working drive
location <- "your_path"
setwd(location)

#Import dataset
dat <- read_xlsx("yourfilepath")

#Explore the data
table(dat$AGEGRP)
table(dat$varstrat)
summary(dat$csp)
summary(dat$HRP2_a)

dat %>% 
  count(PTID) %>% 
  filter(n >1)  # there are no duplicate PTIDs

sum(is.na(dat$ab_wt))

abmiss <- dat %>% 
  filter(is.na(ab_wt))

#Demographics among the children missing sero data
abmiss %>% 
  tabyl(Gender)
abmiss %>% 
  tabyl(AGEGRP) 
abmiss %>% 
  tabyl(AHTYPE)

NG_Ab_data <- dat %>% 
  filter(!is.na(ab_wt)) #For serology analysis, we have 31,234 observations

#Explore the data
skim(NG_Ab_data) #data is very complete
glimpse(NG_Ab_data)

############## Data Cleaning and variable creation #############################

NG_Ab_data$AGEGRP <- as.numeric(NG_Ab_data$AGEGRP) #change AGEGRP from character to numeric variable
NG_Ab_data$AHSTATE <- as.character(NG_Ab_data$AHSTATE)  #change AHSTATE to character for merging

table(NG_Ab_data$A1Q119) 
NG_Ab_data <- NG_Ab_data %>% 
  mutate(net_own= case_when(A1Q119==1~1, #recode net_own variable
                           A1Q119==2~0,
                           A1Q119==8~NA_real_,
                           A1Q119==9~NA_real_),
        Gender= case_when(Gender== 2~ 0, 
                          Gender== 1~ 1), #recode Gender to 0/1
        Urban= case_when(AHTYPE== 2~ 0, 
                         AHTYPE== 1~ 1), #recode AHTYPE to 0/1
        num_nets= ifelse(A1Q119== 2, 0, A1Q120), #recode A1Q120 to have 0 (instead of NA) if HH does not own any nets 
        AGEGRP7= case_when(Age < 3  ~ 1,
                                     3<= Age & Age < 5~ 2,
                                     5<= Age & Age < 7~ 3,
                                     7<= Age & Age < 9~ 4,
                                     9<= Age & Age < 11~ 5,
                                     11<= Age & Age < 13~6,
                                     13<= Age & Age <= 14~7),
        AHLGA2= ifelse(AHLGA < 10, paste0("0", AHLGA), AHLGA),    #add leading 0s before AHLGA < 10
        LGA= str_c(AHSTATE, AHLGA2),   #create unique LGA number
        livestock= case_when(A1Q115==1~1,
                             A1Q115==2~0,
                             A1Q115==8~NA_real_,
                             A1Q115==9~NA_real_),
        aircon= case_when(A1Q105N==1~1,
                        A1Q105N==2~0,
                        A1Q105N==8~NA_real_,
                        A1Q105N==9~NA_real_),
        electric= case_when(A1Q105A==1~1,
                            A1Q105A==2~0,
                            A1Q105A==8~NA_real_,
                            A1Q105A==9~NA_real_),
        water_source= case_when(A1Q101A==1~1,
                                A1Q101A==2~2,
                                A1Q101A==3~0,
                                A1Q101A==8~NA_real_,
                                A1Q101A==9~NA_real_),
        roof_mat= case_when(A1Q108 <= 13~ 0,
                                20< A1Q108  & A1Q108 <23~ 2,
                                30< A1Q108 & A1Q108 <38~1,
                                A1Q108==96~9,
                                A1Q108==97~NA_real_),
        wall_mat= case_when(A1Q109 <= 13~ 0,
                            13< A1Q109  & A1Q109 <26~ 2,
                            30< A1Q109 & A1Q109 <38~1,
                            A1Q109==96~9,
                            A1Q109==97~NA_real_)
         ) 
NG_Ab_data$LGA <- as.numeric(NG_Ab_data$LGA) #convert LGA to numeric

#### Data Exploration ####
NG_Ab_data %>% 
  tabyl(wall_mat, A1Q109)
NG_Ab_data %>% 
  tabyl(A1Q119, net_own)
NG_Ab_data %>% 
  tabyl(net_own)
NG_Ab_data %>% 
  tabyl(Gender) #none missing
NG_Ab_data %>% 
  tabyl(Age) #none missing
NG_Ab_data %>% 
  tabyl(AHTYPE) #none missing
NG_Ab_data %>% 
  tabyl(A1Q119, num_nets)
NG_Ab_data %>% 
  tabyl(AGEGRP)
NG_Ab_data %>% 
  tabyl(PmMSP1pos) #none missing
NG_Ab_data %>% 
  tabyl(PoMSP1pos) #none missing
NG_Ab_data %>% 
  tabyl(PvMSP1pos) #none missing
NG_Ab_data %>% 
  tabyl(HRP2_pos) #112 missing HRP2 antigen positivity 

#Add LGA labels from
source("NGA_labels.r")
NG_Ab_data %>% 
  filter(is.na(LGA)) #none missing an LGA value

distinct(NG_Ab_data, varunit) #3,404 clusters (PSU)
distinct(NG_Ab_data, LGA) #740 LGAs with sero data, 774 LGAs in Nigeria total

NG_Ab_data %>% 
  summarize(n_clusters= n_distinct(varstrat, varunit)) #the clusters are not nested within strata

table(NG_Ab_data$AGEGRP, NG_Ab_data$PmMSP1pos)


my_summary <- NG_Ab_data %>% 
  count(varunit, sort = TRUE)  #number of children per cluster (village)
  my_summary

#Create binary variable for at least 1 mosquito net per 1.8 persons in the Household (0.556 nets/person)
NG_Ab_data <- NG_Ab_data %>% 
  mutate(net_perHHmem = num_nets/ AHMEMBER) %>% 
  mutate(net1.8= ifelse(net_perHHmem >= 0.556, 1, 0))

summary(NG_Ab_data$net_perHHmem)
NG_Ab_data %>% 
  tabyl(net1.8)

#Variable creation for mosquito net coverage
NG_Ab_data <- NG_Ab_data %>% 
  add_count(varunit) %>%  
  rename(n_clust= n)      #n_clust= number of children per cluster
NG_Ab_data <- NG_Ab_data %>% 
  group_by(varunit) %>% 
  add_tally(net1.8) %>% 
  rename(n_net= n) #n_net= number of children per cluster who live in HH with adequate net coverage

#Data Exploration
NG_Ab_data %>% 
  summarize(n_HH= n_distinct(varunit, A1QNUMBER))

new <- NG_Ab_data %>% 
  group_by(varunit) %>% 
  add_count(A1QNUMBER) %>% #number of HHs total= 11,951, and n= number of children per HH
  rename(n_HH= n)    #n_HH= number of children per HH
 
new2 <- new %>% 
  distinct(varunit, A1QNUMBER, .keep_all = TRUE) #dataset of 1 child per HH only
summary(new2$n_HH) #on average, there are 2.6 children per HH
#OR
HH_dataset <- NG_Ab_data %>% distinct(HH_ID, .keep_all = TRUE)
  
new3 <- HH_dataset %>% 
  group_by(varunit) %>% 
  count(HH_ID) %>% 
  add_tally(n) %>% 
  rename(n_clustHH= nn) %>% #n_clustHH= number of HHs per cluster
  distinct(varunit, n_clustHH) 
summary(new3$n_clustHH)

summary(NG_Ab_data$n_clust) 
summary(HH_dataset$net_own) #HH net ownership- 61.8% of HHs own at least one net 


#Create variable for average Mosquito net use in community
NG_Ab_data <- NG_Ab_data %>% 
  mutate(netcover_avg= n_net/ n_clust) %>% 
  mutate(netcover_ind= net1.8- netcover_avg) #added benefit of individual net use compared to community use


#Set options for allowing a single observation per stratum 
options(survey.lonely.psu = "adjust")

#Normalized weights (s_ab_wt)
NG_Ab_data <-  NG_Ab_data %>% 
  group_by(varunit) %>% 
  mutate(weight_mult= n()/ sum(ab_wt)) %>% 
  ungroup() %>% 
  mutate(s_ab_wt= ab_wt* weight_mult)
#Survey Design object
svydat <- svydesign(data=NG_Ab_data, strata=~varstrat, id=~varunit, weights=~s_ab_wt)
summary(svydat)
class(svydat)

#Save dataframe as file
write_xlsx(NG_Ab_data, "Nigeria_Ab_data_Clean.xlsx")

# Mean estimates by State (with St. error and CIs)
svyby(~lsa1, ~AHSTATE, svydat, svymean, vartype=c("se","ci"), na.rm= TRUE)


svymean(~factor(Gender), svydat, na.rm= TRUE) %>%  #overall gender distribution
  confint()
svymean(~factor(AGEGRP), svydat, na.rm= TRUE) %>%  #overall Age group distribution
  confint()
svymean(~factor(AHZONE), svydat, na.rm= TRUE) %>%  #overall Zone distribution
  confint()
svymean(~factor(PmMSP1pos), svydat, na.rm= TRUE) %>%  #Percent positive to PmMSP1
  confint()
svymean(~factor(PoMSP1pos), svydat, na.rm= TRUE) %>%  #Percent positive to PoMSP1
  confint()
svymean(~factor(PvMSP1pos), svydat, na.rm= TRUE) %>%  #Percent positive to PvMSP1
  confint()
svymean(~factor(PfMSP1pos), svydat, na.rm= TRUE) %>%  #Percent positive to PfMSP1
  confint()

##DEMOGRAPHICS N(%) ##
tab1 <- svytable(~Gender, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq)) %>% 
  arrange(desc(Prop)) #arrange highest to lowest
tab1
svytable(~AGEGRP, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq))
svymean(~Age, svydat, na.rm= TRUE)
svyquantile(~Age, svydat, quantiles= c(0,0.25, 0.5, 0.75,1), na.rm= TRUE)

svytable(~Wealthquintile, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq))  #overall Wealth quintile distribution
svytable(~net1.8, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq))  #Net coverage indicator
svytable(~AHZONE, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq))  #Number of children per Zone
svytable(~Urban, svydat) %>% as.data.frame() %>% 
  mutate(Prop= (Freq/ sum(Freq))*100)  #Number of children in urban vs. rural areas
svytable(~AHSTATE, svydat) %>% as.data.frame() %>% 
  mutate(Prop= (Freq/ sum(Freq))*100)  #Number of children per State
svytable(~AHSTATE + PfMSP1pos, svydat) %>% as.data.frame() %>% 
  mutate(Prop= (Freq/ sum(Freq))*100)  #Proportion seropositive to each antigen by State

svytotal(~PmMSP1pos, svydat, na.rm= T) #total weighted num. positive for PmMSP1
svytotal(~PoMSP1pos, svydat, na.rm= T) #total weighted num. positive for PoMSP1
svytotal(~PvMSP1pos, svydat, na.rm= T) #total weighted num. positive for PvMSP1
svytotal(~PfMSP1pos, svydat, na.rm= T) #total weighted num. positive for PfMSP1

#Weighted proportion of children positive by Zone
svyby(~PmMSP1pos, ~AHZONE, svydat, svymean, na.rm= TRUE) %>% 
  confint()
svyby(~PoMSP1pos, ~AHZONE, svydat, svymean, na.rm= TRUE) %>% 
  confint()
svyby(~PvMSP1pos, ~AHZONE, svydat, svymean, na.rm= TRUE) %>% 
  confint()


#### Data prep for basic map of Seroprevalence by state ####

#Weighted proportion of children positive for each antibody by State
pmmsp1_state <- svyby(~PmMSP1pos, ~AHSTATE, svydat, svymean, vartype = ("ci"), na.rm= TRUE)
pomsp1_state <- svyby(~PoMSP1pos, ~AHSTATE, svydat, svymean, vartype = ("ci"), na.rm= TRUE)
pvmsp1_state <- svyby(~PvMSP1pos, ~AHSTATE, svydat, svymean, vartype = ("ci"), na.rm= TRUE)
pfmsp1_state <- svyby(~PfMSP1pos, ~AHSTATE, svydat, svymean, vartype = ("ci"), na.rm= TRUE)


#Save dataframes as Excel files
write_xlsx(pfmsp1_state, "Pf_sero_prev_by_state.xlsx")
write_xlsx(pmmsp1_state, "Pm_sero_prev_by_state.xlsx")
write_xlsx(pomsp1_state, "Po_sero_prev_by_state.xlsx")
write_xlsx(pvmsp1_state, "Pv_sero_prev_by_state.xlsx")

test1 <- merge(pmmsp1_state, pomsp1_state, by= "AHSTATE")
test2 <- merge(pvmsp1_state, pfmsp1_state, by= "AHSTATE")
test3 <- merge(test1, test2, by= "AHSTATE")
merged_data <- test3 %>% 
  select(AHSTATE, PmMSP1pos, PoMSP1pos, PvMSP1pos, PfMSP1pos)
write_xlsx(merged_data, "sero_prev_by_state.xlsx")

#### Data prep for EB Smoothed map of seroprevalence by LGA ####

LGA_n <- svytable(~LGA, svydat) %>% as.data.frame() %>% 
  mutate(Prop= Freq/ sum(Freq))
LGA_n10 <- svytable(~LGA, svydat) %>% as.data.frame() %>% 
  filter(Freq <10)  #54 LGAs have <10 observations

#Weighted proportion of children positive for each antibody by LGA
pmmsp1_LGA <- svyby(~PmMSP1pos, ~LGA, svydat, svymean, vartype = ("se"), na.rm= TRUE) %>% 
    mutate(RSE= (se/ PmMSP1pos)*100)
pomsp1_LGA <- svyby(~PoMSP1pos, ~LGA, svydat, svymean, vartype = ("se"), na.rm= TRUE) %>% 
  mutate(RSE= (se/ PoMSP1pos)*100)
pvmsp1_LGA <- svyby(~PvMSP1pos, ~LGA, svydat, svymean, vartype = ("se"), na.rm= TRUE) %>% 
  mutate(RSE= (se/ PvMSP1pos)*100)
pfmsp1_LGA <- svyby(~PfMSP1pos, ~LGA, svydat, svymean, vartype = ("se"), na.rm= TRUE) %>% 
  mutate(RSE= (se/ PfMSP1pos)*100)
#Number of children per LGA
N_perLGA <- NG_Ab_data %>% 
   tabyl(LGA) 
#Number sero-positive per LGA (weighted)
N_LGA_Pm <- svytable(~LGA+PmMSP1pos, svydat) %>%
  as.data.frame() %>%
  filter(PmMSP1pos== 1)
N_LGA_Po <- svytable(~LGA+PoMSP1pos, svydat) %>%
  as.data.frame() %>%
  filter(PoMSP1pos== 1)
N_LGA_Pv <- svytable(~LGA+PvMSP1pos, svydat) %>%
  as.data.frame() %>%
  filter(PvMSP1pos== 1)
N_LGA_Pf <- svytable(~LGA+PfMSP1pos, svydat) %>%
  as.data.frame() %>%
  filter(PfMSP1pos== 1)

#save dataframes as files
library("openxlsx")
location <- "yourfilepath"
setwd(location)

write.csv(N_perLGA, "N_perLGA.csv")
write.csv(N_LGA_Pf, "N_LGA_pf.csv")
write.csv(N_LGA_Pm, "N_LGA_pm.csv")
write.csv(N_LGA_Po, "N_LGA_po.csv")
write.csv(N_LGA_Pv, "N_LGA_pv.csv")

write.xlsx(pmmsp1_LGA, "Pm_sero_prev_by_LGA.xlsx")
write.csv(pomsp1_LGA, "Po_sero_prev_by_LGA.csv") 
write.csv(pvmsp1_LGA, "Pv_sero_prev_by_LGA.csv") 
write.csv(pfmsp1_LGA, "Pf_sero_prev_by_LGA.csv") 



### MULTI-LEVEL MODELING ####

#Mixed effects, logistic regression exploration
#Nested random effects- Cluster and HH-level
model1a <- glmer(PmMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
              (1 | varunit/ A1QNUMBER),
              data= NG_Ab_data, family= "binomial",
              control=glmerControl(optimizer="bobyqa"))
model1a 
parameters::parameters(model1a, exponentiate = TRUE, details = TRUE)
gof(model1a)

#P. ovale MSP1
model1b <- glmer(PoMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                Urban + (1 | varunit/ A1QNUMBER),
                data= NG_Ab_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa"))
parameters::parameters(model1b, exponentiate = TRUE, details = TRUE)
gof(model1b)
performance::binned_residuals(model1b) #poor fit
#drop Household number Random intercept- only ~2 children per HH on average 

#Final Pm seropositivity model
model1 <- glmer(PmMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                  Urban + (1 | varunit),
                data= NG_Ab_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa"))  
parameters::parameters(model1, exponentiate = TRUE, details = TRUE)
model_performance(model1, metrics= "common")
gof(model1)
performance::check_collinearity(model1)
performance::binned_residuals(model1) 


#Final Po. seropositivity model
model2.2 <- glmer(PoMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                  Urban + (1 | varunit),
                data= NG_Ab_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa"))
parameters(model2.2, exponentiate = TRUE, details = TRUE)
model_performance(model2.2, metrics= "common")
gof(model2.2)
performance::binned_residuals(model2.2) 


#Final Pv seropositivity model
model3 <- glmer(PvMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                  Urban + (1 | varunit),
                data= NG_Ab_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa"))
parameters(model3, exponentiate = TRUE, details = TRUE)
model_performance(model3, metrics= "common") 
gof(model3)
performance::check_collinearity(model3)
performance::binned_residuals(model3) 

#Final Pf seropositivity model
model4 <- glmer(PfMSP1pos ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                  Urban + (1 | varunit/ A1QNUMBER),
                data= NG_Ab_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa"))
parameters(model4, exponentiate = TRUE, details = TRUE)
model_performance(model4, metrics= "common") 
gof(model4)
performance::check_collinearity(model4)



#### SUPPLEMENTARY MATERIALS ####
var_label(NG_Ab_data$HRP2_pos) <- ("HRP2 positive")
val_lab(NG_Ab_data$HRP2_pos)= num_lab( "1 Positive
                                       0 Negative")

# Supplementary Fig. 1- Histogram of children per LGA
LGA_n <- NG_Ab_data %>% 
  group_by(LGA) %>% 
  count(LGA)
summary(LGA_n$n)

h <- hist(LGA_n$n, breaks= 10, xlim=  c(0, 200), xlab= "Number of children per LGA", main= "")

# Supplementary Table 4- Bivariate logistic regression 
summary(NG_Ab_data$AHMEMBER) #on avg. 7 people per house, range: 2- 38

#Repeat for each non-Pf antigen
#Number of HH members
model_Pm1 <- glmer(PmMSP1pos ~ AHMEMBER + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm1, exponentiate = TRUE, details = TRUE)
#livestock
model_Pm2 <- glmer(PmMSP1pos ~ livestock + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm2, exponentiate = TRUE, details = TRUE)

cor.test(NG_Ab_data$Urban, NG_Ab_data$livestock, method= "pearson")
#roof material
model_Pm3 <- glmer(PmMSP1pos ~ factor(roof_mat) + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm3, exponentiate = TRUE, details = TRUE)
#wall material
model_Pm4 <- glmer(PmMSP1pos ~ factor(wall_mat) + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm4, exponentiate = TRUE, details = TRUE)
#Air conditioning
model_Pm5 <- glmer(PmMSP1pos ~ factor(aircon) + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm5, exponentiate = TRUE, details = TRUE)
#Electricity
model_Pm6 <- glmer(PmMSP1pos ~ factor(electric) + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm6, exponentiate = TRUE, details = TRUE)
#Water source
model_Pm7 <- glmer(PmMSP1pos ~ factor(water_source) + (1 | varunit),
                   data= NG_Ab_data, family= "binomial",
                   control=glmerControl(optimizer="bobyqa"))
parameters(model_Pm7, exponentiate = TRUE, details = TRUE)

# Supplementary Fig. 6- Boxplots for each MSP1 assay signal by age
#PmMSP1 boxplot
summary(NG_Ab_data$pmmsp1) #there are some negative numbers for MFI- bg signal, since cannot take the log of a negative #, change values to 1.

NG_Ab_data<- NG_Ab_data %>% 
  mutate(pmmsp1= ifelse(pmmsp1 <= 0, 1, pmmsp1), 
         pomsp1= ifelse(pomsp1 <= 0, 1, pomsp1),
         pvmsp1= ifelse(pvmsp1 <= 0, 1, pvmsp1),
         pfmsp1= ifelse(pfmsp1 <= 0, 1, pfmsp1))

options(scipen= 999)

#Boxplot by all ages- repeat for different antigens
tiff(filename = "Boxplot_PmMSP1_All ages.tiff", width= 2050, height= 1243, units= "px",
     compression = "lzw", res= 300)

pmmsp1_plot <- ggplot(NG_Ab_data, aes(x= factor(Age), y=pmmsp1)) + 
  geom_boxplot()
pmmsp1_plot +  
  scale_y_log10()+ 
  labs(x= "Age (years)", y= "Assay signal (MFI- bg)", title= "Assay signal for PmMSP1")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black")) 
dev.off()


# Supplementary Fig. 8- Air conditioning and urban/rural LOESS curves
data <- NG_Ab_data %>% 
  filter(Age >0 & Urban== 1 & aircon== 1) #Update for each combination of Urban/aircon

#Get dataframes with seroprevalence for each age group for plotting points
#Air conditioning and Urban
df_PmMSP17 <- data %>% 
  tabyl(AGEGRP7, PmMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0")
colnames(df_PmMSP17) <- c("age", "PmMSP1pos") 
#Need to plot AGEGRP at the mid-point of each age group
df_new_Pm7 <-edit(df_PmMSP17) 
print(df_new_Pm7)

df_PoMSP17 <- data %>% 
  tabyl(AGEGRP7, PoMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0")
colnames(df_PoMSP17) <- c("age", "PoMSP1pos") 
#Need to plot AGEGRP at the mid-point of each age group
df_new_Po7 <-edit(df_PoMSP17) 
print(df_new_Po7)

df_PvMSP17 <- data %>% 
  tabyl(AGEGRP7, PvMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0")
colnames(df_PvMSP17) <- c("age", "PvMSP1pos") 
#Need to plot AGEGRP at the mid-point of each age group
df_new_Pv7 <-edit(df_PvMSP17) 
print(df_new_Pv7)

#Seroprevalence by all ages for LOESS regression
df_PmSP1 <- data %>% 
  tabyl(Age, PmMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0") #take out the 0 column
colnames(df_PmMSP1) <- c("age", "PmMSP1pos") 

df_PoMSP1 <- data %>% 
  tabyl(Age, PoMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0") #take out the 0 column
colnames(df_PoMSP1) <- c("age", "PoMSP1pos") 

df_PvMSP1 <- data %>% 
  tabyl(Age, PvMSP1pos) %>% 
  adorn_percentages(denominator = "row") %>% 
  dplyr::select(-"0") #take out the 0 column
colnames(df_PvMSP1) <- c("age", "PvMSP1pos") 

#LOESS
loessMod50 <- loess(PvMSP1pos ~ age, data=df_PvMSP1, span=0.50) # 50% smoothing span
smoothed50 <- predict(loessMod50, se= TRUE) 

#Plotting function
plot_loess <- function(X, Y, title, color){
  plot(X, Y, main= title, 
       xlab="Age", ylab="Proportion Seropositive",
       pch=16, cex=1.1,
       xlim=c(0,14), ylim=c(0,0.2), col= color)
}

ages = seq(from = 0, to=14,by=0.1)
preds <- predict(loessMod50, newdata = data.frame(age= ages), 
                 interval = 'confidence')

grid.newpage()
par(mfrow=c(1,1))
#print the figure
tiff(filename = "LOESS_PvMSP1_All ages.tiff", width= 1700, height= 879, units= "px",
     compression = "lzw", res= 200)

plot_loess(X= df_new_Pv7$age, Y= df_new_Pv7$PvMSP1pos, title= "PvMSP1", "blue")
lines(smoothed50$fit, x= df_PvMSP17$age, col="blue") #Urban + A/C

points(x= df_new_Pv7$age, y= df_new_Pv7$PvMSP1pos,  pch= 16, cex=1.1, col="skyblue")
lines(smoothed50$fit, x= df_PvMSP1$age, col="skyblue")  #Urban, no A/C

points(x= df_new_Pv7$age, y= df_new_Pv7$PvMSP1pos,  pch= 16, cex=1.1, col="red")
lines(smoothed50$fit, x= df_PvMSP1$age, col="red")  #Rural + A/C

points(x= df_new_Pv7$age, y= df_new_Pv7$PvMSP1pos,  pch= 16, cex=1.1, col="lightcoral")
lines(smoothed50$fit, x= df_PvMSP1$age, col="lightcoral") #Rural, no A/C

legend(x= "topleft", lty= 1, lwd= 2, col= c("blue", "skyblue", "red", "lightcoral"), text.col= "black",
       legend= c("Urban & A/C", "Urban, no A/C", "Rural & A/C", "Rural, no A/C"))
dev.off()

### Supplementary Fig. 9- Boxplot MSP1 signal by PF HRP2 antigen Positivity
AB_antigen_data <- NG_Ab_data %>% 
  filter(!(is.na(HRP2_pos)))
library(patchwork)
library(rstatix)
library(ggpubr)
library("rcompanion")

plot_pm1 <- ggplot(AB_antigen_data, aes(x= factor(HRP2_pos), y=pmmsp1)) + 
  geom_boxplot() +  
  scale_y_log10()+ 
  labs(x= "P. falciparum infection", y= "Assay signal (MFI- bg)", title= "Assay signal for PmMSP1")+
  theme(plot.title = element_text(face= "bold", size= 10, hjust= 0.5, color= "black")) 

plot_po1 <- ggplot(AB_antigen_data, aes(x= factor(HRP2_pos), y=pomsp1)) + 
  geom_boxplot() +  
  scale_y_log10()+ 
  labs(x= "P. falciparum infection", y= "Assay signal (MFI- bg)", title= "Assay signal for PoMSP1")+
  theme(plot.title = element_text(face= "bold", size= 10, hjust= 0.5, color= "black")) 

plot_pv1 <- ggplot(AB_antigen_data, aes(x= factor(HRP2_pos), y=pvmsp1)) + 
  geom_boxplot()+
  scale_y_log10()+ 
  labs(x= "P. falciparum infection", y= "Assay signal (MFI- bg)", title= "Assay signal for PvMSP1")+
  theme(plot.title = element_text(face= "bold", size= 10, hjust= 0.5, color= "black")) 
#print figure
tiff(filename = "Boxplot_MSP1 v Pfinfection.tiff", width= 1400, height= 777, units= "px",
     compression = "lzw", res= 140)

plot_pm1 + plot_po1 + plot_pv1
dev.off()



