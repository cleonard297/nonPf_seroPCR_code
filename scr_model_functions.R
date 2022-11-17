
#####################################################################
#####################################################################
##                                                                 ##
## CODE FOR FITTING CATALYTIC CONVERSION MODELS TO  ANTIBODY DATA  ##
##                                                                 ##
##  Colleen Leonard                                                ##
##  1/12/2022                                                      ##
##                                                                 ##
## Code from: Michael White                                        ##
## Imperial College London                                         ##
## m.white08@imperial.ac.uk                                        ##
##  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6303365/          ##
##  Model 1, *assumes constant transmission                        ##
#####################################################################
#####################################################################

## call the libraries that are needed for analysis
library(tidyverse)
library(grid)
library(gridExtra)
library(MASS)
library(compiler)
library(binom)

###0.1 Data preparation- Done in scr_models.R and scr_models_supp_materials.R ###


###############################################
## 0.2 bin the data by age group to plot
## define the bins for your data depending
## on the age ranges available in your data

age_bins     <- seq(from=0, to=14, by=2)				  ## age bins will be set from the youngest age group to the maximum
age_bins_mid <- seq(from=1.5, to=13.5, by=2) 		  
age_bins

N_bins <- length(age_bins) - 1 #Number of bins= 7 
N_bins

SP_range_bins <- matrix(NA, nrow=N_bins, ncol=3)  
colnames(SP_range_bins) <- c("seroprev", "low", "high")	  ## generate an empty matrix to be filled with sero-prev for each age bin and the binomial CIs 

#For getting the binomial confidence intervals
for(i in 1:N_bins){
  
  index <- intersect( which(AB_data[,1]>age_bins[i]), which(AB_data[,1]<=age_bins[i+1]) ) #inclusive- first group= [0-10]
  temp  <- AB_data[index,2]
  
  SP_range_bins[i,] <- as.numeric(as.vector(
    binom.confint( sum(AB_data[index,2]), length(index), method="wilson")[1,4:6]
  ))
}


############## 1.0 MODEL    ############################## 
## 1.1 specify the serological model which may
## have generated our sero-prev data 


par_MC <- c(0.1, 0.5)  ## (lambda, rho)

model_M1 <- function(a, par)
{
  lambda <- par[1]
  rho    <- par[2]
  
  SP_prop <- ( lambda/(lambda+rho) )*( 1 - exp(-(lambda+rho)*a) )
  
  SP_prop
}

model_M1 <- cmpfun(model_M1, options=list(optimize=3)) 

###################################################
## 1.2 calculate the likelihood of the paramater set
## using a binomial likelihood
## specify the fixed value of the SRR

loglike_M1 <- function(par)
{
  lambda <- par[1]
  rho    <- par[2]
  
  SP_model <- model_M1( AB_data[,1], par )
  
  mu <- log(SP_model) 
  
  loglike <- AB_data[,2]*log(SP_model) + (1-AB_data[,2])*log(1-SP_model)
  
  sum( loglike )
}

loglike_M1 <- cmpfun(loglike_M1, options=list(optimize=3)) ## use the loglike function above to calculate the likelihhod for a specific parameter set


###################################################
## 1.3 Define the priors for each estimated parameter
## we currently define uniform priors for all parameters

LARGE = 1e10     			## Large value for rejecting parameters with prior

prior_M1 <- function( par )
{
  lambda <- par[1]
  rho    <- par[2]
  
  ######################################
  ## Uniform prior on lambda ~ U(0,10)
  
  if( lambda>0 && lambda<10 )
  {
    prior_lambda <- log(1/10)
  }else{
    prior_lambda <- -LARGE
  }
  
  ######################################
  ## Uniform prior on rho ~ U(0,10)
  
  if( rho>0 && rho<10 )
  {
    prior_rho <- log(1/10)
  }else{
    prior_rho <- -LARGE
  }
  
  prior <- prior_lambda + prior_rho
  
  prior
}

prior_M1 <- cmpfun(prior_M1, options=list(optimize=3)) ## generate the priors needed for each MCMC evaluation 

############## 2.0 MCMC    ##############################



N_mcmc       <- 10000      ## Number of MCMC iterations
N_tune_start <- 300        ## Start of adaptive tuning of covariance matrix of MVN proposals
N_tune_end   <- 3000       ## End of adaptive tuning of covariance matrix of MVN proposals
N_adapt      <- 4000       ## End of adaptive scaling of proposal size with rm_scale in 2.1

step_scale  <- 1           ## Scaler for step size
MCMC_accept <- 0           ## Track the MCMC acceptance rate

max_corr    <- 0.75        ## Maximum degree of correlation


#################################################
## 2.1 Robbins-munro step scaler


rm_scale <- function(step_scale, mc, log_prob)
{
  dd <- exp(log_prob)
  if( dd < -30 ){ dd <- 0 }
  dd <- min( dd, 1 )
  
  rm_temp <- ( dd - 0.23 )/( (mc+1)/(0.01*N_adapt+1) )
  
  out <- step_scale*exp(rm_temp)
  
  out <- max( out, 0.05 )
  out <- min( out, 5)
  out
}


#################################################
## 2.2 Prepare object to store MCMC fitting output

MCMC_par           <- matrix(NA, nrow=N_mcmc, ncol=4)
colnames(MCMC_par) <- c("lambda", "rho", "loglike", "prior")


#########################################################
## 2.3 Implement MCMC iterations

par_MC <- c(0.1, 0.5)  ## (lambda, rho)	 ## initial paramater guess

Sigma_MC <- diag( (0.25*par_MC)^2 )      ## Initial guess of covariance of MVN proposal dist



loglike_MC <- loglike_M1( par_MC ) + prior_M1( par_MC )



for(mc in 1:N_mcmc)
{
  par_MCp1 <- mvrnorm(n=1, mu=par_MC, Sigma=step_scale*Sigma_MC)
  
  
  if( prior_M1(par_MCp1) > -0.5*LARGE  ){
    
    loglike_MCp1 <- loglike_M1( par_MCp1 ) + prior_M1( par_MCp1 )
    
    
    log_prob <- min( loglike_MCp1-loglike_MC, 0 )           
    
    if( log(runif(1)) < log_prob ) 
    {
      par_MC <- par_MCp1
      
      loglike_MC  <- loglike_MCp1
      MCMC_accept <- MCMC_accept + 1                       
    }
    
    #######################################
    ## RM scaling of proposal step size
    
    if( mc < N_adapt ){
      step_scale <- rm_scale( step_scale, mc, log_prob)
    }
    
    #######################################
    ## Adaptive tuning of covariance matrix
    
    if( (mc > N_tune_start) && (mc < N_tune_end) )
    {
      cov_MC <- cov( MCMC_par[1:(mc-1),1:2] )
      
      ###########################
      ## Checks for tuning
      
      if( min(diag(cov_MC)) > 1e-6 )
      {
        ###########################
        ## Check for high degree of correlation
        
        sd_MC_inv <- 1/sqrt(diag(cov_MC))
        
        corr_MC <- t(t(cov_MC*sd_MC_inv)*sd_MC_inv)
        
        corr_MC[intersect( which( corr_MC > max_corr ), which(corr_MC<0.99999) )] <- max_corr
        corr_MC[  which( corr_MC < -max_corr )] <- -max_corr
        
        t(t(corr_MC*(1/sd_MC_inv))*(1/sd_MC_inv))
        
        Sigma_MC <- cov_MC
      }
    }
  }
  
  MCMC_par[mc,1:2] <- par_MC				## store the output parameter values from each MCMC realisation 
  MCMC_par[mc,3]   <- loglike_MC			## return the associated loglikelihood
  MCMC_par[mc,4]   <- prior_M1(par_MC)	## return the prior
}



#########################################################
## 2.4 Examine MCMC chains

par(mfrow=c(1,3))



#####################################
## PANEL 1: alpha MCMC chain Trace plot

 plot(x=1:N_mcmc, y=MCMC_par[,1], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="lambda", 
     main="lambda")


#####################################
## PANEL 2: rr MCMC chain Trace Plot

plot(x=1:N_mcmc, y=MCMC_par[,2], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="rho", 
     main="rho")

#####################################
## PANEL 3: likelihood Trace Plot

plot(x=1:N_mcmc, y=MCMC_par[,3], 
     pch=19, col="grey", cex=0.25,
     xlab="MCMC iteration", ylab="likelihood", 
     ylim=quantile(MCMC_par[,3], prob=c(0.01,1)),
     main="likelihood" )


#########################################################
## 2.5 Examine posterior distribution

MCMC_burn <- MCMC_par[floor(0.2*nrow(MCMC_par)):(nrow(MCMC_par)-1),]

par(mfrow=c(1,2))

#####################################
## PANEL 1: alpha MCMC posterior (SCR) Density plot


DEN_SCR <- density( MCMC_burn[,1] )

QUANT_SCR <- quantile( MCMC_burn[,1], prob=c(0.025, 0.5, 0.975) )


plot(x=DEN_SCR$x, y=DEN_SCR$y, type='l',
     xlim=c(0, max(DEN_SCR$x)),
     xlab="lambda", ylab="", 
     main="posterior: lambda" )


low_index  = which(DEN_SCR$x<QUANT_SCR[1])
mid_index  = intersect( which(DEN_SCR$x>=QUANT_SCR[1]), which(DEN_SCR$x<=QUANT_SCR[3]) )
high_index = which(DEN_SCR$x>QUANT_SCR[3])

polygon( x=c( DEN_SCR$x[low_index], rev(DEN_SCR$x[low_index]) ),
         y=c( rep(0,length(low_index)), rev(DEN_SCR$y[low_index]) ), 
         col="pink")

polygon( x=c( DEN_SCR$x[mid_index], rev(DEN_SCR$x[mid_index]) ),
         y=c( rep(0,length(mid_index)), rev(DEN_SCR$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN_SCR$x[high_index], rev(DEN_SCR$x[high_index]) ),
         y=c( rep(0,length(high_index)), rev(DEN_SCR$y[high_index]) ), 
         col="pink")

points(x=rep(QUANT_SCR[2],2), y=c(0,max(DEN_SCR$y)), type='l', lty="dashed", lwd=2)





#####################################
## PANEL 2: rr MCMC posterior (SRR) Density plot


DEN_RR <- density(MCMC_burn[,2] )

QUANT_RR <- quantile( MCMC_burn[,2], prob=c(0.025, 0.5, 0.975) )

plot(x=DEN_RR$x, y=DEN_RR$y, type='l',
     xlim=c(0, max(DEN_RR$x)),
     xlab="rho", ylab="", 
     main="posterior: rho" )


low_index  = which(DEN_RR$x<QUANT_RR[1])
mid_index  = intersect( which(DEN_RR$x>=QUANT_RR[1]), which(DEN_RR$x<=QUANT_RR[3]) )
high_index = which(DEN_RR$x>QUANT_RR[3])

polygon( x=c( DEN_RR$x[low_index], rev(DEN_RR$x[low_index]) ),
         y=c( rep(0,length(low_index)), rev(DEN_RR$y[low_index]) ), 
         col="pink")

polygon( x=c( DEN_RR$x[mid_index], rev(DEN_RR$x[mid_index]) ),
         y=c( rep(0,length(mid_index)), rev(DEN_RR$y[mid_index]) ), 
         col="grey")

polygon( x=c( DEN_RR$x[high_index], rev(DEN_RR$x[high_index]) ),
         y=c( rep(0,length(high_index)), rev(DEN_RR$y[high_index]) ), 
         col="pink")

points(x=rep(QUANT_RR[2],2), y=c(0,max(DEN_RR$y)), type='l', lty="dashed", lwd=2)


################ 3.0 TEST ###################

#############################################
## 3.1 Extract posterior medians and 
##     calculate model prediction


par_median <- apply(X=MCMC_burn[,1:2], MARGIN=2, FUN = median)
  par_lower_ci  <- apply(X=MCMC_burn[,1:2], MARGIN=2, FUN = quantile, prob= 0.025) 
  par_upper_ci <- apply(X=MCMC_burn[,1:2], MARGIN=2, FUN = quantile, prob= 0.975)

age_seq <- seq(from=0, to=14, by=1)

M1_predict <- model_M1(age_seq, par_median )
M1_lowerl <- model_M1(age_seq, par_lower_ci)
M1_upperl <- model_M1(age_seq, par_upper_ci)

