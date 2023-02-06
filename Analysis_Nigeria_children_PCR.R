#******************************************************************#
#      ANALYSIS OF NIGERIA CHILDRENS MALARIA PET-PCR DATA          #
#******************************************************************#
# author: Colleen Leonard
# date: January 25, 2022

 
#Load required packages
library(haven)
library(plyr) 
library(rgdal)
library(lme4)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
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
library(detectseparation) #checking a glmm for complete separation
library(numDeriv)
library(see)
library(car) #identifying outliers and problems with residuals
library(brglm) #Firth's penalized-likelihood logistic regression for small sample bias

# set working drive
location <- "yourworkingdirectory"
setwd(location)

#Import raw Ab dataset, including all 32,494 samples (225 samples with GST > 500)
dat <- read_xlsx("Nigeria_children_sero_Total.xlsx")
#Import PCR dataset
PCR_dat <- read_xlsx("Nigeria-Merge with PET data_clean.xlsx")

NG_Ab_data <- dat %>% 
  mutate(net_own= case_when(A1Q119==1~1, #recode net_own variable
                            A1Q119==2~0,
                            A1Q119==8~NA_real_,
                            A1Q119==9~NA_real_),
         Urban= case_when(AHTYPE== 2~ 0, 
                          AHTYPE== 1~ 1), #recode AHTYPE to 0/1
         num_nets= ifelse(A1Q119== 2, 0, A1Q120))   #Create binary variable for at least 1 mosquito net per 1.8 persons in the Household (0.556 nets/person)

NG_Ab_data <- NG_Ab_data %>% 
  mutate(net_perHHmem = num_nets/ AHMEMBER) %>% 
  mutate(net1.8= ifelse(net_perHHmem >= 0.556, 1, 0))
NG_Ab_data <- NG_Ab_data %>% 
  add_count(varunit) %>%  
  rename(n_clust= n)      #n_clust= number of children per cluster
NG_Ab_data <- NG_Ab_data %>% 
  group_by(varunit) %>% 
  add_tally(net1.8) %>% 
  rename(n_net= n) #n_net= number of children per cluster who live in HH with adequate net coverage

#Create variable for average Mosquito net coverage community- Created from the total sample of 32,269
NG_Ab_data <- NG_Ab_data %>% 
  mutate(netcover_avg= n_net/ n_clust) %>% 
  mutate(netcover_ind= net1.8- netcover_avg) #added benefit of individual net use compared to community use

#Limit Sero data to only variables needed to add to PCR dataset
NG_Ab_data_filt <- NG_Ab_data %>% 
  select(PTID, HH_ID, varunit, A1QNUMBER, AHTYPE, Urban, Wealthquintile, netcover_avg, netcover_ind)

#Check if the PTID columns in PCR_dat are the same- they are. 
identical(PCR_dat$PTID...1, PCR_dat$PTID...5)

PCR_data= subset(PCR_dat, select= -c(PTID...5)) #drop the repeat PTID var

PCR_data= rename(PCR_data, PTID= PTID...1)  #rename variable 

#Filter the PCR dataset a bit
PCR_data <- PCR_data %>% 
  subset(select= -c(agegp, Guspec, guspec, Plate_num_ab, Order_ab, Order_ag, Well_ab, Well_ag)) %>% 
  filter(is.na(`missing samples`)) #if sample was missing remove from dataset

#Join the PCR data with the filtered Ab Data
NG_PCR_data <- PCR_data %>% 
  left_join(NG_Ab_data_filt, "PTID")

#Create AGEGRP category for all samples
NG_PCR_data <- NG_PCR_data %>% 
  mutate(AGEGRP= case_when(Age < 5~ 1,
                           5<= Age & Age < 10~ 2,
                           10<= Age & Age < 15~ 3))
table(NG_PCR_data$Age, NG_PCR_data$AGEGRP)

skim(NG_PCR_data) #No covariates missing from the Full (Ab) dataset
NG_PCR_data <- NG_PCR_data %>% 
  mutate(Gender= case_when(Gender== 2~ 0, 
                           Gender== 1~ 1)) #recode Gender to 0/1

NG_PCR_data %>% 
    get_dupes(HH_ID)  #about 200 households with at least 2 children. 

NG_PCR_data %>% 
  count(varunit) %>%   
  filter(n>4)  #only 5 clusters with 5+ people

identical(NG_PCR_data$cluster, NG_PCR_data$varunit) #cluster number from PCR data is different from the varunit (cluster var from Sero data.)

test <- NG_PCR_data %>% 
  group_by(cluster) %>% 
  add_count(cluster) %>% 
  rename(n_clust= n) %>%  #number of children per cluster
  select(cluster, n_clust, A1QNUMBER, HH_ID)
#if children are in the same "varunit" they are also in the same "cluster" although the numbers are different.
clust <- test %>% 
  distinct(cluster, .keep_all= TRUE)
summary(clust$n_clust) #on average, there are 1.4 persons per Cluster

#Create new variables for PCR positive to Pm, Po and Pv
NG_PCR_data <- NG_PCR_data %>% 
  mutate(PCR_pf= case_when(`species present`==1 | 
                          `species present` =="2 (Pf Pm)" |
                            `species present`== "2 (Pf PO)" |
                            `species present`== "2(Pf Pm)" |
                            `species present`== "3 (Pf Pm PO)"|
                            `species present`== "3 (Pf PM PO)" ~1, 
                            `species present`== "1 (Pm)"|
                            `species present`== "1 (PO)"|
                            `species present`== "PO"|
                            `species present`== "Pm"|
                            `species present`== "2 (Pm PO)" ~0),
         
         PCR_pm= case_when( `species present`== "1 (Pm)"|
                               `species present`== "Pm"|
                               `species present` =="2 (Pf Pm)" |
                               `species present`== "2(Pf Pm)" |
                               `species present`== "3 (Pf Pm PO)"|
                               `species present`== "3 (Pf PM PO)" | 
                               `species present`== "2 (Pm PO)" ~1,
                               `species present`== 1 |
                               `species present`== "1 (PO)"|
                               `species present`== "PO"|
                               `species present`== "2 (Pf PO)" ~0),
         PCR_po= case_when(  `species present`== "1 (PO)"|
                                    `species present`== "PO"|
                                    `species present`== "3 (Pf Pm PO)"|
                                    `species present`== "3 (Pf PM PO)" | 
                                    `species present`== "2 (Pm PO)"|
                                    `species present`== "2 (Pf PO)" ~1,
                                    `species present`== 1 |
                                    `species present`== "Pm"|
                                    `species present` =="2 (Pf Pm)" |
                                    `species present`== "2(Pf Pm)" |
                                    `species present`== "1 (Pm)"~0))

NG_PCR_data %>% 
  tabyl(`species present`, PCR_pf)

NG_PCR_data <- mutate_at(NG_PCR_data, c("PCR_pf", "PCR_pm", "PCR_po"), ~replace(., is.na(.), 0))   #NAs are negative by PCR, so should be 0


NG_PCR_data <- NG_PCR_data %>% 
  mutate(PCR_species= ifelse(is.na(`species present`), 0, `species present`))

NG_PCR_data %>% 
  tabyl(`species present`, PCR_species)

NG_PCR_data <-  NG_PCR_data %>% 
  mutate(PCR_pfmono= ifelse(PCR_species== "1", 1, 0))
NG_PCR_data %>% 
  tabyl(PCR_species, PCR_pfmono)

NG_PCR_data <-  NG_PCR_data %>% 
  mutate(PCR_pm_pf= ifelse(PCR_species== "1 (Pm)" |
                             PCR_species== "Pm" | 
                             PCR_species== "2 (Pf Pm)"|
                             PCR_species== "2(Pf Pm)", 1, 0),
         PCR_po_pf= ifelse(PCR_species== "2 (Pf PO)" |
                             PCR_species== "1 (PO)"|
                             PCR_species== "PO", 1, 0))
NG_PCR_data %>% 
  tabyl(PCR_species, PCR_pm_pf)
NG_PCR_data %>% 
  tabyl(PCR_species, PCR_po_pf)

#Exporting Complete PCR Dataset
write_xlsx(NG_PCR_data, "Nigeria_PCR_data_complete.xlsx")

#Data exploration
NG_PCR_data %>% 
  tabyl(AGEGRP, PCR_pm) %>% 
  adorn_percentages(denominator = "row") %>% 
  adorn_pct_formatting(digits= 1)

NG_PCR_data %>% 
  tabyl(AGEGRP, PCR_po) %>% 
  adorn_percentages(denominator = "row") %>% 
  adorn_pct_formatting(digits= 1)

#### Multi-level Models- model fitting ####
#P. malariae Models- postive for Pm. by PCR vs. not 
model1a <- glmer(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                   Urban + (1 | cluster/ A1QNUMBER),
                data= NG_PCR_data, family= "binomial",
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
parameters::parameters(model1a, exponentiate = TRUE, details = TRUE)
summary(model1a)
#Checking for complete separation 
fixed_form= PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind + Urban
glm(fixed_form, data = NG_PCR_data, family = binomial, method="detect_separation") 
#There is no complete separation for Pm

NG_PCR_data %>% 
  tabyl(Wealthquintile, PCR_pm) #only 9 people with Pm. in the wealthiest category
xtabs(~PCR_pm + AGEGRP, data= NG_PCR_data)
xtabs(~PCR_pm + Gender, data= NG_PCR_data)
xtabs(~PCR_pm + Urban, data= NG_PCR_data) #all represented by an adequate number of people
NG_PCR_data %>%
  subset(PCR_pm==1) %>%
  tabyl(Wealthquintile, Urban)
NG_PCR_data %>%
  subset(PCR_po==1) %>%
  tabyl(Wealthquintile, Urban)
#For the outcome PCR Pm= yes, there is only 1 person where Wealth=5 AND Urban= 0
#For the outcome PCR Po= yes, there are no persons in the Wealthiest category 
 
fixed_form= PCR_po ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind + Urban
glm(fixed_form, data = NG_PCR_data, family = binomial, method="detect_separation") 
#Confirmed there is complete separation for Po Wealth variable

#Since there are only on average 1.3 persons per Cluster, median= 1 person per cluster- can ignore clustering 
#Source: https://jech.bmj.com/content/62/8/752.long

model2 <- glm(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                Urban,
              data= NG_PCR_data, family= "binomial")
summary(model2) #the deviance residuals appear to be not symmetrical
performance::binned_residuals(model2) 
parameters::parameters(model2, exponentiate = TRUE, details = TRUE)
model_performance(model2, metrics= "common")
gof(model2)  # the model fits the data well
performance_hosmer(model2) #***Model does not fit well
performance::check_collinearity(model2) #no multicollinearity

#Try using Firth's penalized likelihood with wealth categories, since there's only one person in wealthiest/rural category

#P. malariae Model- Better model fit
model_pm3 <- brglm(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                    netcover_avg, data= NG_PCR_data, family= "binomial", method= "brglm.fit")
performance::binned_residuals(model_pm3) 
parameters::parameters(model_pm3, exponentiate = TRUE, details = TRUE)
model_performance(model_pm3, metrics= "common")
gof(model_pm3)
performance_hosmer(model_pm3) #***Model seems to fit well by Hosmer-Lemeshow GOF Test

#P. ovale Models 
model4 <- glm(PCR_po ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + netcover_avg + netcover_ind +
                Urban,
              data= NG_PCR_data, family= "binomial")
summary(model4) 
parameters::parameters(model4, exponentiate = TRUE, details = TRUE)

performance::binned_residuals(model4) #Only 60% of the residuals are inside the error bounds, probably bad model fit
model_performance(model4, metrics= "common")
gof(model4)  
performance::check_collinearity(model4)

residualPlots(model4) #none of the variables show a pattern of their residuals 
outlierTest(model4) #No outliers as determined by Bonferroni p
#data exploration
NG_PCR_data %>% 
  tabyl(PCR_po)  #Only 42 children with P. ovale! 
NG_PCR_data %>% 
  tabyl(Wealthquintile, PCR_po) #0 people with Po. in the wealthiest category
xtabs(~PCR_po + AGEGRP, data= NG_PCR_data)
xtabs(~PCR_po + Gender, data= NG_PCR_data)
xtabs(~PCR_po + Urban, data= NG_PCR_data) #all represented by an adequate number of people

# Since we have only ~40 events for P. ovale, we have a small sample  bias
##Use firth's penalized likelihood to reduce the bias

#P. ovale Model- better fit 
model_po_1 <- brglm(PCR_po ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                      netcover_avg, data= NG_PCR_data, family= "binomial", method= "brglm.fit")
summary(model_po_1)
parameters::parameters(model_po_1, exponentiate = TRUE, details = TRUE)
performance::binned_residuals(model_po_1)
gof(model_po_1) #model fits the data
performance_hosmer(model_po_1)


# Supplementary table 2.0- Outcome: Pm. (and Po.) infection, using Pf mono-infection as the Reference group

#Create a temporary dataset which only includes Pf, and any Pm infections (mono or mixed)
PmPf_dat <- NG_PCR_data %>% 
  filter(PCR_pm== 1 | PCR_pfmono== 1)
PmPf_dat %>% 
  tabyl(PCR_pm)
model_pm2.0 <- brglm(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                     netcover_avg, data= PmPf_dat, family= "binomial", method= "brglm.fit")

performance::binned_residuals(model_pm2.0) #97% of the residuals are inside the error bounds
parameters::parameters(model_pm2.0, exponentiate = TRUE, details = TRUE)
model_performance(model_pm2.0, metrics= "common")
gof(model_pm2.0)
performance_hosmer(model_pm2.0)

#Create a temporary dataset which only includes Pf mono-infection, and Any Po infection
PoPf_dat <- NG_PCR_data %>% 
  filter(PCR_po== 1 | PCR_pfmono== 1)

PoPf_dat %>% 
  tabyl(PCR_species, PCR_po)

model_po2.0 <- brglm(PCR_po ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= PoPf_dat, family= "binomial", method= "brglm.fit")

performance::binned_residuals(model_po2.0)
parameters::parameters(model_po2.0, exponentiate = TRUE, details = TRUE)
model_performance(model_po2.0, metrics= "common")
gof(model_po2.0)
performance_hosmer(model_po2.0)


#3.0 Using PCR negative as the Reference Group

#Create a temporary dataset which only includes Pm (mixed or mono) and PCR negative samples
Po_neg_dat <- NG_PCR_data %>% 
  filter(PCR_po== 1 | PCR_species== 0)

Pm_neg_dat %>% 
  tabyl(PCR_pm, Wealthquintile) #at least 9 in each category

model_pm3.0 <- brglm(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Pm_neg_dat, family= "binomial", method= "brglm.fit")

performance::binned_residuals(model_pm3.0) #100% of the residuals are inside the error bounds
parameters::parameters(model_pm3.0, exponentiate = TRUE, details = TRUE)
model_performance(model_pm3.0, metrics= "common")
performance_hosmer(model_pm3.0)

#Create a temporary dataset which only includes Po (mixed or mono) and PCR negative samples
Po_neg_dat <- NG_PCR_data %>% 
  filter(PCR_po== 1 | PCR_species== 0)

Po_neg_dat %>% 
  tabyl(PCR_po, Wealthquintile) #no obs. with Po in wealthcat= 5

model_po3.0 <- brglm(PCR_po ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Po_neg_dat, family= "binomial", method= "brglm.fit")

performance::binned_residuals(model_po3.0)
parameters::parameters(model_po3.0, exponentiate = TRUE, details = TRUE)
model_performance(model_po3.0, metrics= "common")
performance_hosmer(model_po3.0)


#Create a temporary dataset which only includes Pf (mixed or mono) and PCR negative samples
Pf_neg_dat <- NG_PCR_data %>% 
  filter(PCR_pf== 1 | PCR_species== 0)

model_pf3.0 <- brglm(PCR_pf ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Pf_neg_dat, family= "binomial", method= "brglm.fit")

performance::binned_residuals(model_pf3.0)
parameters::parameters(model_pf3.0, exponentiate = TRUE, details = TRUE)
model_performance(model_pf3.0, metrics= "common")
performance_hosmer(model_pf3.0)


#### 4.0 FINAL MODELS- Table 2- Any Pm (or Po) Infection vs. Negative (by all 4 antigens) ####
#Read in cleaned Antigen dataset (NG_antigen_data)

#Create a temporary dataset which only includes Pm (mixed or mono) and Antigen negative samples
#PCR abbreviated data set
NG_PCR_abrv <- NG_PCR_data %>% 
  select(PTID, PCR_species, PCR_pm, PCR_pf, PCR_po, PCR_pm_pf, PCR_po_pf, PCR_pfmono)

#Check if all PCR obs (1204) are found within the full antigen dataset= Yes
check <- do.call(paste0, NG_PCR_abrv[,1]) %in% do.call(paste0, NG_antigen_data[,1]) 
check[1000:1204]

#Merge all Antigen with PCR data
NG_PCRant_dat <- NG_antigen_data %>% 
  left_join(NG_PCR_abrv, "PTID")

#Filter for Pm PCR positive
PCR_pm_dat <- NG_PCRant_dat %>% 
  filter(PCR_pm==1)
Neg_4ant <- NG_PCRant_dat %>% 
  filter(HRP2_pos== 0, pLDH_pos== 0, PvLDH_pos== 0, pAldo_pos== 0)

Pm_PCRant_dat <- rbind(Neg_4ant, PCR_pm_dat) %>% 
  arrange(desc(PCR_pm)) %>% 
  distinct(PTID, .keep_all = TRUE)  #there are some duplicate obs (9) which were PCR_pm positive, but negative for all 4 antigens 
 #Keeping obs where PCR_pm= 1

Pm_PCRant_dat <- mutate_at(Pm_PCRant_dat, c("PCR_pm"), ~replace(., is.na(.), 0))   #NAs are negative by PCR

write_xlsx(Pm_PCRant_dat, "Nigeria_Pm_PCR_ant_dat.xlsx")

model_pm4.0 <- brglm(PCR_pm ~ factor(AGEGRP) + Gender + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Pm_PCRant_dat, family= "binomial", method= "brglm.fit")

parameters::parameters(model_pm4.0, exponentiate = TRUE, details = TRUE)
performance::binned_residuals(model_pm4.0)
model_performance(model_pm4.0, metrics= "common")
gof(model_pm4.0)
performance_hosmer(model_pm4.0) #model seems to fit well


#Filter for Po PCR positive
PCR_po_dat <- NG_PCRant_dat %>% 
  filter(PCR_po==1)

Po_PCRant_dat <- rbind(Neg_4ant, PCR_po_dat) %>% 
  arrange(desc(PCR_po)) %>% 
  distinct(PTID, .keep_all = TRUE)

Po_PCRant_dat <- mutate_at(Po_PCRant_dat, c("PCR_po"), ~replace(., is.na(.), 0))   #NAs are negative by PCR

write_xlsx(Po_PCRant_dat, "Nigeria_PO_PCR_ant_dat.xlsx")

Po_PCRant_dat %>% 
  tabyl(Wealthquintile, PCR_po)

model_po4.0 <- brglm(PCR_po ~ Gender + factor(AGEGRP) + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Po_PCRant_dat, family= "binomial", method= "brglm.fit")

parameters::parameters(model_po4.0, exponentiate = TRUE, details = TRUE)
performance::binned_residuals(model_po4.0)
model_performance(model_po4.0, metrics= "common")
gof(model_po4.0)
performance_hosmer(model_po4.0) #model seems to fit well


#Filter for Pf PCR positive
PCR_pf_dat <- NG_PCRant_dat %>% 
  filter(PCR_pf==1)

Pf_PCRant_dat <- rbind(Neg_4ant, PCR_pf_dat) %>% 
  arrange(desc(PCR_pf)) %>% 
  distinct(PTID, .keep_all = TRUE)

Pf_PCRant_dat <- mutate_at(Pf_PCRant_dat, c("PCR_pf"), ~replace(., is.na(.), 0))   #NAs are negative by PCR
Pf_PCRant_dat %>% 
  tabyl(PCR_pf) #over 5% positive for Pf 
Pf_PCRant_dat %>% 
  tabyl(PCR_pf, Wealthquintile)

write_xlsx(Pf_PCRant_dat, "Nigeria_Pf_PCR_ant_dat.xlsx")

#regular GLM 
model_pf4.0 <- glm(PCR_pf ~ Gender + factor(AGEGRP) + factor(Wealthquintile) + Urban + netcover_ind +
                       netcover_avg, data= Pf_PCRant_dat, family= "binomial")

parameters::parameters(model_pf4.0, exponentiate = TRUE, details = TRUE)
performance::binned_residuals(model_pf4.0)
model_performance(model_pf4.0, metrics= "common")
gof(model_pf4.0)
performance_hosmer(model_pf4.0) #model seems to fit well



#### Map of PCR data ####
#Read in cleaned PCR dataset
NG_PCR_data <- read_xlsx("PCR_data_complete.xlsx")

#Number of PCR samples per State, raw
N_perState <- NG_PCR_data %>% 
  tabyl(AHSTATE)

PM_state <- NG_PCR_data %>% 
  tabyl(AHSTATE, PCR_pm)
colnames(PM_state) <- c("AHSTATE", "PM_neg", "PM_pos") 

Po_by_state <- NG_PCR_data %>% 
  tabyl(AHSTATE, PCR_po)
colnames(Po_by_state) <- c("AHSTATE", "PO_neg", "PO_pos") 

Pf_by_state <- NG_PCR_data %>% 
  tabyl(AHSTATE, PCR_pf)
colnames(Pf_by_state) <- c("AHSTATE", "PF_neg", "PF_pos") 

#merge dataframes
test1 <- merge(PM_state, Po_by_state, by= "AHSTATE")
PCR_by_State <- merge(test1, Pf_by_state, by= "AHSTATE")
PCR_by_State <- PCR_by_State %>% 
  mutate(total_n= PM_neg + PM_pos)
#save PCR data by state

#PCR samples per State, by selection strategy
NG_PCR_data %>% 
  tabyl(`Selection Criteria`)
NG_PCR_data <- NG_PCR_data %>% 
  mutate(selection_cri= case_when(`Selection Criteria`== "Zones"~ 1,
                                  `Selection Criteria`== "pAldolase" |
                                  `Selection Criteria`=="pLDH" |
                                  `Selection Criteria`== "pLDH+/HRP2-" |
                                   `Selection Criteria`== "PvLDH"~2,
                                  ))

Selection <- NG_PCR_data %>% 
  tabyl(AHSTATE,selection_cri)
write.csv(Selection, "Zones_Sel_State.csv")

PM_sel_State <- NG_PCR_data %>% 
  tabyl(AHSTATE,selection_cri, PCR_pm)
Pm_pos <- PM_sel_State[2] #positives only
write_xlsx(Pm_pos, "Pm_pos.xlsx")

Po_sel_State <- NG_PCR_data %>% 
  tabyl(AHSTATE,selection_cri, PCR_po)
Po_pos <- Po_sel_State[2]
write_xlsx(Po_pos, "Po_pos.xlsx")

Pf_sel_State <- NG_PCR_data %>% 
  tabyl(AHSTATE,selection_cri, PCR_pf)
Pf_pos <- Pf_sel_State[2]
write_xlsx(Pf_pos, "Pf_pos.xlsx")


#### SUPPLEMENTARY MATERIALS ####
#Variable label
var_label(NG_PCR_data$AGEGRP) <- ("Age Group")
#Value labels
val_lab(NG_PCR_data$PCR_pm)= num_lab( "1 Positive
                                       0 Negative")
val_lab(NG_PCR_data$PCR_po)= num_lab( "1 Positive
                                       0 Negative")
val_lab(NG_PCR_data$AGEGRP)= num_lab( "1 0- 4 yrs
                                       2 5- 9 yrs    
                                       3 10- 14 yrs")


### Supplementary Fig. 9- Boxplot MSP1 signal by PF HRP2 antigen Positivity
NG_PCR_data %>% 
  filter(pmmsp1 < 0.5)
NG_PCR_data<- NG_PCR_data %>% 
  mutate(pmmsp1= ifelse(pmmsp1 <= 0, 1, pmmsp1), 
         pomsp1= ifelse(pomsp1 <= 0, 1, pomsp1),
         pvmsp1= ifelse(pvmsp1 <= 0, 1, pvmsp1),
         pfmsp1= ifelse(pfmsp1 <= 0, 1, pfmsp1))


options(scipen= 999)


### Supplementary Fig. 10
#Save Boxplot- PmMSP1 vs. Pm. PCR 
tiff(filename = "PmMSP1 by Pm PCR_Age boxplots.tiff", width= 1650, height= 820, units= "px",
     compression = "lzw", res= 200)

plot_pm1 <- ggplot(NG_PCR_data, aes(x= factor(PCR_pm), y=pmmsp1, color= factor(AGEGRP))) + 
  geom_boxplot()
plot_pm1 +  
  scale_y_log10()+ 
  labs(x= "P. malariae mono or mixed infection (PCR)", y= "Assay signal (MFI- bg)", title= "Assay signal for PmMSP1", color= "Age Group")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black"))+
  scale_fill_discrete(name= "Age Group")

dev.off()

#Save Boxplot- PmMSP1 vs. Po. PCR 
tiff(filename = "PmMSP1 by Po PCR_Age boxplots.tiff", width= 1650, height= 820, units= "px",
     compression = "lzw", res= 200)
plot_pm12 <- ggplot(NG_PCR_data, aes(x= factor(PCR_po), y=pmmsp1, color= factor(AGEGRP))) + 
  geom_boxplot()
plot_pm12 +  
  scale_y_log10()+ 
  labs(x= "P. ovale mono or mixed infection (PCR)", y= "Assay signal (MFI- bg)", title= "Assay signal for PmMSP1", color= "Age Group")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black")) 

dev.off()

#Save Boxplot- PoMSP1 vs. Po. PCR 
tiff(filename = "PoMSP1 by Po PCR_Age boxplots.tiff", width= 1650, height= 820, units= "px",
     compression = "lzw", res= 200)
plot_pm2 <- ggplot(NG_PCR_data, aes(x= factor(PCR_po), y=pomsp1, color= factor(AGEGRP))) + 
  geom_boxplot()
plot_pm2 +  
  scale_y_log10()+ 
  labs(x= "P. ovale mono or mixed infection (PCR)", y= "Assay signal (MFI- bg)", title= "Assay signal for PoMSP1", color= "Age Group")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black")) 
dev.off()

#Save Boxplot- PoMSP1 vs. Pm. PCR 
tiff(filename = "PoMSP1 by Pm PCR_Age boxplots.tiff", width= 1650, height= 820, units= "px",
     compression = "lzw", res= 200)
plot_pm2.1 <- ggplot(NG_PCR_data, aes(x= factor(PCR_pm), y=pomsp1, color= factor(AGEGRP))) + 
  geom_boxplot()
plot_pm2.1 +  
  scale_y_log10()+ 
  labs(x= "P. malariae mono or mixed infection (PCR)", y= "Assay signal (MFI- bg)", title= "Assay signal for PoMSP1", color= "Age Group")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black")) 
dev.off()


### Supplementary Fig. 11- Scatter plot 
SCR_bystate <- read_xlsx("SCR_by_State.xlsx")

#lm equation
install.packages("ggpubr")
library(ggpubr)

tiff(filename = "Scatterplot_Pm SCR vs inf prev.tiff", width= 2000, height= 1230, units= "px",
     compression = "lzw", res= 300)
plot1 <- ggplot(SCR_bystate, aes(x= Pm_PCR_prev, y=SCR_pm)) + 
  geom_point()
plot1 +  
  geom_smooth(method = lm)+
  stat_regline_equation(label.y= 0.195, aes(label = c(..rr.label.., ..p.value.label..)))+
  labs(x= "P. malariae mono or mixed infection prevalence (%)", y= "Annual rate of seroconversion")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black"))+
  theme_classic()

dev.off()

#P. ovale Scatterplot
tiff(filename = "Scatterplot_Po SCR vs inf prev.tiff", width= 1900, height= 1048, units= "px",
     compression = "lzw", res= 300)
plot2 <- ggplot(SCR_bystate, aes(x= Po_PCR_prev, y=SCR_po)) + 
  geom_point()
plot2 +  
  xlim(0, 6)+
  geom_smooth(method = lm)+
  stat_regline_equation(label.y= 0.056, aes(label = ..rr.label..))+
  labs(x= "P. ovale mono or mixed infection prevalence (%)", y= "Annual rate of seroconversion")+
  theme(plot.title = element_text(face= "bold", size= 12, hjust= 0.5, color= "black"))+
  theme_classic()

dev.off()


