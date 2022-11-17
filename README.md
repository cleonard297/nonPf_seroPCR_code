# nonPf_seroPCR_code

# Background: 
 This repository houses the code used to analyze non-falciparum malaria serology and PCR data collected during the 2018 nationwide Nigeria HIV/AIDS Indicator and Impact Survey (NAIIS) household survey. Demographic, socioeconomic, and behavioral risk factors were explored in relation to malaria infection and exposure using multilevel logistic regression. Seroconversion models were used to model the rate of change from seronegative to seropositive for each of the three non-falciparum antigens (PmMSP1, PoMSP1, and PvMSP1).  

Code files:
Data cleaning and logistic regression: 
  	 NAIIS_data_cleaning.sas
 	 Analysis_Nigeria_children_sero.R
   Analysis_Nigeria_children_PCR.R 
   Analysis_Nigeria_children_alldata.R

Empirical bayes smoothed maps:
 	EB_smoothing.R

Seroconversion models:
  Scr_model_functions.R
  Scr_models.R
  Scr_models_by_State.R
