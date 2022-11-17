####Code for assessing SCR curves by covariates ####
#read in Antibody data
NG_Ab_data <- read_xlsx("yourdatafile.xlsx")
data <- NG_Ab_data %>% 
  filter(Age >0 & Gender== 1)  #Change covariate combinations of interest

data %>% 
  tabyl(PvMSP1pos, AGEGRP7) %>% 
  adorn_percentages(denominator ="col")
data %>% 
  tabyl(Age, AGEGRP7)

Kin_pre_binned <- cbind(data$Age, data$PvMSP1pos)		  ## make an array of all individuals age and sero-status for analysis 

colnames(Kin_pre_binned) <- c("age", "PvMSP1pos_bin") 

AB_data <- Kin_pre_binned  							   ## assign the object to a new name that will be evaluated by the function  


location <- "yourworkingdirectory"
setwd(location)

source("scr_model_functions.r")

QUANT_SCR  #Seroconversion estimate
QUANT_RR   #Seroreversion estimate
library(coda)

#Re-run the scr_model_functions for the additional covariate categories 

###############################################
## Plot data and model prediction SCR curves

# Supplementary Figure 7- SCR for each antigen by Gender
grid.newpage()
par(mfrow=c(1,1))
#Change scale depending on Species
pmscale <- 0.8
poscale <- 0.4
pvscale <- 0.25

#tiff("PvMSP1 SCR models By Gender.tiff", units= "px", width= 1225, height= 650, res= 300)
plot(x=age_bins_mid, y=SP_range_bins[,1], 
     pch=16, cex=1.1,
     xlim=c(0,14), ylim=c(0,0.2), col= "red",
     xlab="Age (years)", ylab="Proportion seropositive", 
     main="P. vivax MSP1" )

for(i in 1:N_bins) #For CI error bars
{
  arrows(x0=age_bins_mid[i], y0=SP_range_bins[i,2], 
         x1=age_bins_mid[i], y1=SP_range_bins[i,3], 
         length=0.03, angle=90, code=3, col= "red", lwd=1)	
}

points(x=age_seq, y=M1_predict, 
       type='l', lwd=2, col="red")  

polygon(c(age_seq,rev(age_seq)),c(M1_lowerl,rev(M1_upperl)),col=adjustcolor("red",alpha.f = 0.25), border=NA)

#Function for Next Lines
graph <- function(scale, color){
points(x=age_bins_mid, y=SP_range_bins[,1], pch=16, cex=1.1, col= color,
       xlim=c(0,14), ylim=c(0,scale))
for(i in 1:N_bins) #For CI error bars
{
  arrows(x0=age_bins_mid[i], y0=SP_range_bins[i,2], 
         x1=age_bins_mid[i], y1=SP_range_bins[i,3], 
         length=0.03, angle=90, code=3, col=color, lwd=1)	
}
points(x=age_seq, y=M1_predict, 
       type='l', lwd=2, col=color)  
polygon(c(age_seq,rev(age_seq)),c(M1_lowerl,rev(M1_upperl)),col=adjustcolor(color,alpha.f = 0.25), border=NA)
}

#NEXT line
graph(scale=0.2, color= "blue")  #Male

legend(x= "topright", lty= 1, lwd= 2, col= c("blue", "red"), text.col= "black",
      legend= c("Male", "Female"))

#dev.off() #removes all plots and Saves the high res tiff


# Code for Supplementary Figure 11- SCR for each antigen by State

data <- NG_Ab_data %>% 
  filter(Age >0 & AHSTATE=="35") #Change number for each State

Kin_pre_binned <- cbind(data$Age, data$PoMSP1pos)		  ## make an array of all individuals age and sero-status for analysis 

colnames(Kin_pre_binned) <- c("age", "Po_MSP1_bin") 

AB_data <- Kin_pre_binned  							   ## assign the object to a new name that will be evaluated by the function  

source("scr_model_functions.r")

QUANT_SCR  #Seroconversion estimate
QUANT_RR   #Seroreversion estimate
#Re-run the scr_model_functions for the additional States




