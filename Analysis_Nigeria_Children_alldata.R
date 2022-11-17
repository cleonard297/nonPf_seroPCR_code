#Load packages
library(dplyr)
library(readxl)
library(writexl)
library(janitor)
library(tidyverse)
library(expss)
library(survey)


#Read in cleaned Antigen dataset 
NG_antigen_data <- read_xlsx("Nigeria_children_Antigen_abbrev.xlsx")
NG_antigen_data %>% 
  tabyl(HRP2_pos)

#Preparing full antigen dataset for logistic regression
############## Data Cleaning and variable creation #############################
NG_antigen_data$AGEGRP <- as.numeric(NG_antigen_data$AGEGRP) #change AGEGRP from character to numeric variable
NG_antigen_data$AHSTATE <- as.character(NG_antigen_data$AHSTATE)  #change AHSTATE to character for merging

table(NG_antigen_data$A1Q119) #only 5 cases where mosquito net information is missing.
NG_antigen_data <- NG_antigen_data %>% 
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
         LGA= str_c(AHSTATE, AHLGA2)
  ) 
NG_antigen_data$LGA <- as.numeric(NG_antigen_data$LGA) #convert LGA to numeric

NG_antigen_data %>% 
  tabyl(A1Q119, net_own)
NG_antigen_data %>% 
  tabyl(net_own)
NG_antigen_data %>% 
  tabyl(Gender) #none missing
NG_antigen_data %>% 
  tabyl(Age) #none missing


distinct(NG_antigen_data, varunit) #3,407 clusters (PSU) with Antigen data
distinct(NG_antigen_data, LGA) #741 LGAs with Antigen data, 774 LGAs in Nigeria total


#Create binary variable for at least 1 mosquito net per 1.8 persons in the Household (0.556 nets/person)
NG_antigen_data <- NG_antigen_data %>% 
  mutate(net_perHHmem = num_nets/ AHMEMBER) %>% 
  mutate(net1.8= ifelse(net_perHHmem >= 0.556, 1, 0))
NG_antigen_data %>% 
  tabyl(net1.8)

NG_antigen_data <- NG_antigen_data %>% 
  add_count(varunit) %>%  
  rename(n_clust= n)      #n_clust= number of children per cluster
NG_antigen_data <- NG_antigen_data %>% 
  group_by(varunit) %>% 
  add_tally(net1.8) %>% 
  rename(n_net= n) #n_net= number of children per cluster who live in HH with adequate net coverage

#Create variable for average Mosquito net use in community
NG_antigen_data <- NG_antigen_data %>% 
  mutate(netcover_avg= n_net/ n_clust) %>% 
  mutate(netcover_ind= net1.8- netcover_avg)

#Set options for allowing a single observation per stratum 
options(survey.lonely.psu = "adjust")

#Normalized weights (s_ag_wt)
NG_antigen_data <-  NG_antigen_data %>% 
  group_by(varunit) %>% 
  mutate(weight_mult= n()/ sum(ag_wt)) %>% 
  ungroup() %>% 
  mutate(s_ag_wt= ag_wt* weight_mult)
#Survey Design object
svydat <- svydesign(data=NG_antigen_data, strata=~varstrat, id=~varunit, weights=~s_ag_wt)
class(svydat)


#### Number of samples with Antigen data per State (2nd selection strategy for PCR assay) ####

N_per_state <- NG_antigen_data %>% 
  tabyl(AHSTATE)
getwd()
write_xlsx(N_per_state, "N_per_State_antigen.xlsx")

#Number HRP2 positive per state
PF_by_state <-NG_antigen_data %>% 
  tabyl(AHSTATE, HRP2_pos)
colnames(PF_by_state) <- c("States", "Pf Neg", "Pf Positive")
getwd()
write_xlsx(PF_by_state, "Pf_per_State_HRP2.xlsx")


