######################################################################################################################
# Author: Alex Keil
# Program: g_computation.R
# Language: R
# Date: August 2022
# Project: Environmental justice bootcamp
# Description: Carry out g-computation analysis of a plasmode simulated dataset
#  to demonstrate a method for estimating impacts of exposure inequity. See
#  accompanying slides for details
#
# User notes: the only necessary change you will need to make is on line 21, 
#  where you will need to specify the path to the directory in which 
#  plasmode_data_nhanes.csv is stored
#
# Keywords: environmental justice, causal inference, disparities
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

library(readr)
library(dplyr)
library(ggplot2)

# CHANGE ME TO A LOCAL PATH WHERE THE DATA ARE STORED
pth = "~/temp/ej/data/"
setwd(pth)


# select analytic variables (ignore NHANES weights for now)
dat <- read_csv("plasmode_data_nhanes.csv", show_col_types = FALSE)

##############################################################################
# Basic notes:
#  g-computation involves 2 primary steps: 
#    1) fitting models for the data
#    2) making predictions from those models
#  Many details go into those steps, which are given below. Sometimes I leave
#   in a few different options for accomplishing similar tasks in hopes that
#   this code can be adapted for different uses.
#
#  Please send questions or comments to alex.keil@nih.gov
# 
##############################################################################

##############################################################################
# Programming notes:
# Much of the computation occurs within custom made functions written below
#  This approach can hurt readability, but is useful for code useability
#  and adaptability to new problems. The alternative to just code everything
#  directly without functions is possible, but is more error prone and is
#  harder to adapt because every variable name would need to be changed to 
#  use on new data.
##############################################################################

##############################################################################
##############################################################################
# step 1: fit a model for the outcome
##############################################################################
##############################################################################

# r-style formula (linear function of exposures)
#  See nhanes descriptions of each of these variable names (2013-2014 data)
# https://wwwn.cdc.gov/Nchs/Nhanes/search/datapage.aspx?Component=Examination&CycleBeginYear=2013
# "lbxhscrp_plasmode" is a simulated outcome based on lbxhscrp
# "female" is obtained by simply subtracting 1 from the gender variable
# "dmdeduc2_alt" groups value 7 (uknown) with value 5 (college or above)
#   This was done for convenience rather than principle (but it was a single observation)
#   and can be considered a single imputation approach
form1 = lbxhscrp_plasmode ~ urxop1 + urxop2 + lbxbcd + lbxbpb +
  ridageyr + I(ridageyr^2) + 
  bmxbmi + I(bmxbmi^2) + lbxcot +
  factor(ridreth1)+
  factor(dmdeduc2_alt)+
  female
fit1 <- glm(form1, data=dat, family = gaussian())
summary(fit1)
##############################################################################
##############################################################################
# step 2: create "intervention" datasets to reflect the contrasts of interest
##############################################################################
##############################################################################

# natural course (no intervention) - possibly using large Monte Carlo sample (see note below in step 3)
#dat_nc = dat[sample(1:nrow(dat), size=50000, replace=TRUE),]
#dat_nc = do.call(rbind, replicate(100, dat, simplify=FALSE)) # alternative to reduce randomness: make many copies of the observed data
dat_nc = dat # for a linear model, this is actually all that is needed (which makes this code very fast)
# data in which Black non-Hispanic participants have same geometric mean as White non-Hispanic participants
##############################################################################
# Section notes:
# intervention method used here: "stochastic intervention"
#   post-intervention exposure based off of proportional change in the observed exposure
#  where the proportion is based on geometric mean ratios. An "intervention" 
#  in g-computation (for the time-fixed case) is a simple as creating a new
#  dataset with manipulated exposure values. Here, 4 exposures are manipulated in 
#  the dataset "dat_wnh_stochastic"
##############################################################################


# getting geometric means
geometric_mean = function(x,...) exp(mean(log(x),...))
# joint "intervention" on both exposures of interest: need to 
#  first identify geometric mean exposures for racial/ethnic groups of interest
geomean_wnh1 = with(dat[dat$ridreth1==3,], geometric_mean(urxop1))
geomean_bnh1 = with(dat[dat$ridreth1==4,], geometric_mean(urxop1))
geomean_wnh2 = with(dat[dat$ridreth1==3,], geometric_mean(urxop2))
geomean_bnh2 = with(dat[dat$ridreth1==4,], geometric_mean(urxop2))
geomean_wnhcd = with(dat[dat$ridreth1==3,], geometric_mean(lbxbcd))
geomean_bnhcd = with(dat[dat$ridreth1==4,], geometric_mean(lbxbcd))
geomean_wnhpb = with(dat[dat$ridreth1==3,], geometric_mean(lbxbpb))
geomean_bnhpb = with(dat[dat$ridreth1==4,], geometric_mean(lbxbpb))

dat_wnh_stochastic = dat_nc %>%
  mutate(
    urxop1 = case_when(
      # new exposure for Black NH participants is obtained by dividing the original exposure by the geometric mean ratio
      # improve identifiability by ensuring that intervention exposures are in the same support as the observed exposures (i.e. they have the same LOD)
      ridreth1 == 4 ~ pmax(urxop1 / (geomean_bnh1/geomean_wnh1), min(urxop1)),
      TRUE ~ urxop1
    ),
    urxop2 = case_when(
      ridreth1 == 4 ~ pmax(urxop2 / (geomean_bnh2/geomean_wnh2), min(urxop2)),
      TRUE ~ urxop2
    ),
    lbxbcd = case_when(
      ridreth1 == 4 ~ pmax(lbxbcd / (geomean_bnhcd/geomean_wnhcd), min(lbxbcd)),
      TRUE ~ lbxbcd
    ),
    lbxbpb = case_when(
      ridreth1 == 4 ~ pmax(lbxbpb / (geomean_bnhpb/geomean_wnhpb), min(lbxbpb)),
      TRUE ~ lbxbpb
    )
  )
# confirm geometric means are (nearly) the same across Black/White participants
geometric_mean(dat_wnh_stochastic[dat_wnh_stochastic$ridreth1==3,]$urxop1)
geometric_mean(dat_wnh_stochastic[dat_wnh_stochastic$ridreth1==4,]$urxop1)




##############################################################################
##############################################################################
# step 3: predict outcomes from original model using intervention data sets
##############################################################################
##############################################################################


##############################################################################
# Section notes: 
#  generally, we need Monte Carlo sample to get mean outcomes under intervention.
#  Here, that is not needed because the outcome is untransformed and a linear 
#  model is used, but code to carry out the MC sampling is given in the previous
#  section
#  
##############################################################################
set.seed(1232) # for replicability only: results in practical examples should be robust to the seed value
#errsd =  predict(fit1, se.fit=TRUE)$residual.scale
errsd =  0  # this is a programming trick to make it easy to adapt to using a Monte Carlo
            # sample if needed. Just comment out this line and uncomment the previous line
ncpreds =   (predict(fit1, newdat= dat_nc) +             rnorm(nrow(dat_nc),0, errsd))
int1preds = (predict(fit1, newdat= dat_wnh_stochastic) + rnorm(nrow(dat_nc),0,errsd))
# natural course mean should be similar to orignal mean (may differ some according to simulation error)
#mean(ncpreds)
#mean(dat$lbxhscrp_plasmode)

##############################################################################
##############################################################################
# step 4: effect measure point estimates
##############################################################################
##############################################################################

##############################################################################
# Section notes: 
#  This uses the model predictions to calculate a population averaged additive 
#  change in CRP upon modifying the exposure distribution of Black, 
#  non-Hispanic participants to marginally resemble the exposure distribution 
#  of the White non-Hispanic participants
##############################################################################
bnhidx <- which(dat_nc$ridreth1 == 4)
wnhidx <- which(dat_nc$ridreth1 == 3)
# effect among Black, non-Hispanic participants
ncmean_bnh <- mean(ncpreds[bnhidx])
ncmean_wnh <- mean(ncpreds[wnhidx])
meandiff_bnh <- mean(int1preds[bnhidx]) - mean(ncpreds[bnhidx])
# effect in the population after only intervening among Black, non-Hispanic participants
meandiff_marginal <- mean(int1preds) - mean(ncpreds)
(nc_crp_disparity <- ncmean_bnh - ncmean_wnh)
# this will be identical to meandiff_bnh when there is no intervention among the
#  White NH participants
diff_crp_disparity <- (mean(int1preds[bnhidx]) - mean(int1preds[wnhidx])) - (nc_crp_disparity)

  
##############################################################################
##############################################################################
# step 4: non-parametric bootstrap for confidence interval estimation
##############################################################################
##############################################################################

##############################################################################
# Section notes: 
#  This just repeats everything that was done above, but does so repeatedly.
#  An efficient (and less readible) way to write code for g-computation is
#  to start out with something like the "bootiter" function below, and
#  just apply it to your original data instead of bootstrap resamples.
#  As a side effect, that approach is less prone to bugs (so if you do
#  see differences between the bootstrap approach and steps 1-3 above,
#  it may be a bug)
##############################################################################

# a function to perform a g-computation analysis of a single bootstrap dataset
bootiter <- function (dat, modelform=form1){
  # resample with replacement
  bootdat <- dat[sample(1:nrow(dat), replace=TRUE),]
  # if using MC sampling, use one of the two following lines instead of just setting dat_nc_boot = bootdat
  #dat_nc_boot = bootdat[sample(1:nrow(bootdat), size=50000, replace=TRUE),]
  #dat_nc_boot = do.call(rbind, replicate(100, bootdat, simplify=FALSE))
  dat_nc_boot = bootdat
  # fit model
  fitboot <- glm(modelform, data=bootdat, family = gaussian())
  # joint "intervention" on both exposures of interest
  # re-estimate geometric means (implies slightly different intervention,
  #  which accounts for the uncertainty in the intervention based on 
  #  sampling variability)
  geomean_wnh1 =  with(bootdat[bootdat$ridreth1==3,], geometric_mean(urxop1))
  geomean_bnh1 =  with(bootdat[bootdat$ridreth1==4,], geometric_mean(urxop1))
  geomean_wnh2 =  with(bootdat[bootdat$ridreth1==3,], geometric_mean(urxop2))
  geomean_bnh2 =  with(bootdat[bootdat$ridreth1==4,], geometric_mean(urxop2))
  geomean_wnhcd = with(bootdat[bootdat$ridreth1==3,], geometric_mean(lbxbcd))
  geomean_bnhcd = with(bootdat[bootdat$ridreth1==4,], geometric_mean(lbxbcd))
  geomean_wnhpb = with(bootdat[bootdat$ridreth1==3,], geometric_mean(lbxbpb))
  geomean_bnhpb = with(bootdat[bootdat$ridreth1==4,], geometric_mean(lbxbpb))
  
  dat_wnh_stochastic_boot = dat_nc_boot %>%
    mutate(
      urxop1 = case_when(
        ridreth1 == 4 ~ pmax(urxop1 / (geomean_bnh1/geomean_wnh1), min(urxop1)),
        TRUE ~ urxop1
      ),
      urxop2 = case_when(
        ridreth1 == 4 ~ pmax(urxop2 / (geomean_bnh2/geomean_wnh2), min(urxop2)),
        TRUE ~ urxop2
      ),
      lbxbcd = case_when(
        ridreth1 == 4 ~ pmax(lbxbcd / (geomean_bnhcd/geomean_wnhcd), min(lbxbcd)),
        TRUE ~ lbxbcd
      ),
      lbxbpb = case_when(
        ridreth1 == 4 ~ pmax(lbxbpb / (geomean_bnhpb/geomean_wnhpb), min(lbxbpb)),
        TRUE ~ lbxbpb
      )
    )
  #ncpreds_boot = predict(fit1, newdat= bootdat)
  #int1preds_boot = predict(fit1, newdat= dat_wnh_stochastic_boot)
  #errboot = predict(fitboot, se.fit=TRUE)$residual.scale
  errboot = 0
  ncpreds_boot = (predict(fitboot, newdat= dat_nc_boot) + rnorm(nrow(dat_nc_boot),0, errboot))
  int1preds_boot = (predict(fitboot, newdat= dat_wnh_stochastic_boot) + rnorm(nrow(dat_nc_boot),0, errboot))
  #
  bnhidx <- which(dat_nc_boot$ridreth1 == 4)
  wnhidx <- which(dat_nc_boot$ridreth1 == 3)
  # effect among Black, non-Hispanic participants
  ncmean_bnh <- mean(ncpreds_boot[bnhidx])
  ncmean_wnh <- mean(ncpreds_boot[wnhidx])
  meandiff_bnh <- mean(int1preds_boot[bnhidx]) - mean(ncpreds_boot[bnhidx])
  # effect in the population after only intervening among Black, non-Hispanic participants
  meandiff_marginal <- mean(int1preds_boot) - mean(ncpreds_boot)
  nc_crp_disparity <- mean(ncpreds_boot[bnhidx]) - mean(ncpreds_boot[wnhidx])
  diff_crp_disparity <- (mean(int1preds_boot[bnhidx]) - mean(int1preds_boot[wnhidx])) - (nc_crp_disparity)
  c(ncmean_wnh=ncmean_wnh, 
    ncmean_bnh=ncmean_bnh, 
    meandiff_bnh=meandiff_bnh,
    meandiff_marginal=meandiff_marginal, 
    nc_crp_disparity=nc_crp_disparity,
    diff_crp_disparity = diff_crp_disparity
    )
}

# a wrapper function that just runs the bootstrap analysis on a bunch of datasets
# and compiles the results into a matrix
dobootstrap <- function(dat, numboots=2000){
  resList = list()
  for(iter in 1:numboots){
    if((iter %% 100)==0) cat(".")
    resList[[iter]] = bootiter(dat, modelform=form1)
  }
  cat("\n")
  do.call(rbind, resList)
}


# now run the bootstrap in a loop (there are often more efficient ways, but this is clearer)
set.seed(1232)
boot_estimates = dobootstrap(dat, numboots=2000)
# a matrix of bootstrap results: column 1 is "meandiff_bnh" and column 2 is "meandiff_marginal"


# a few check on the bootstrap iterations before believing them
# check for bias in bootstrap results, if so, use bias corrected bootstrap
(point_estimates <- c(ncmean_wnh, ncmean_bnh, meandiff_bnh, meandiff_marginal, nc_crp_disparity, diff_crp_disparity)) # point estimates
(bootmeans <- apply(boot_estimates, 2, mean))            # bootstrap means
(bootbias <- bootmeans-point_estimates)

# check for skewness in bootstrap results, if so, consider using bias-corrected, accelerated bootstrap (can be tricky to implement, so it's out of scope here)
hist(boot_estimates[,1], breaks=30, main=colnames(boot_estimates)[1])
hist(boot_estimates[,2], breaks=30, main=colnames(boot_estimates)[2])
hist(boot_estimates[,3], breaks=30, main=colnames(boot_estimates)[3])
hist(boot_estimates[,4], breaks=30, main=colnames(boot_estimates)[4])
hist(boot_estimates[,5], breaks=30, main=colnames(boot_estimates)[5])
hist(boot_estimates[,6], breaks=30, main=colnames(boot_estimates)[6])

# a function with 3 methods for calculating bootstrap confidence intervals:
#  percentile based, bias-corrected percentile based, and normality based
#  Typically, bias-corrected percentile based will be most general and useful
bootstrap_ci <- function(pointestimates, bootestimates, alpha=0.05){
  nvars = length(pointestimates)
  pbootlessthanE <- numeric(nvars)
  for(i in 1:nvars){
    # bias on the cumulative distribution function basis (assuming normality)
    pbootlessthanE[i] = mean(bootestimates[,i] < pointestimates[i])
  }
  W = qnorm(pbootlessthanE, mean = 0, sd = 1)
  oldq = c(alpha/2, 1-alpha/2)
  z_old = qnorm(oldq)
  newq = cbind(pnorm(z_old[1] + W*2),pnorm(z_old[2] + W*2)) # quantiles of bootstrap distribution to use to yield unbiased intervals
  percentile_intervals = matrix(nrow=nvars, ncol=2)
  bc_intervals = matrix(nrow=nvars, ncol=2)
  norm_intervals = matrix(nrow=nvars, ncol=2)
  dimnames(percentile_intervals)[2] = list(c("LCI_pctl", "UCI_pctl"))
  dimnames(bc_intervals)[2] = list(c("LCI_bcpctl", "UCI_bcpctl"))
  dimnames(norm_intervals)[2] = list(c("LCI_normal", "UCI_normal"))
  stderr = apply(bootestimates, 2, sd)
  for(i in 1:nvars){
    percentile_intervals[i,] = quantile(bootestimates[,i], oldq)
    bc_intervals[i,] = quantile(bootestimates[,i], newq[i,])
    norm_intervals[i,] = pointestimates[i] + z_old* stderr[i]
  }
  list(Std.Err=stderr,bc_intervals=bc_intervals, percentile_intervals=percentile_intervals, norm_intervals=norm_intervals)
}


bootci_estimates = bootstrap_ci(pointestimates = point_estimates, bootestimates=boot_estimates)

# a function to compile bootstrap results with original estimates and print them nicely
printEffectEstimates <- function(pointestimates, bootciestimates){
  se = bootciestimates$Std.Err
  pval = 2-2*pnorm(abs(pointestimates/se))
  printCoefmat(cbind(Estimate=pointestimates, Std.Err=se, 
                     do.call(cbind, bootciestimates[-c(1)]), `Pr(>|z|)`=pval), 
               signif.stars = FALSE)
}

##############################################################################
# print final estimates
##############################################################################

printEffectEstimates(
  pointestimates = point_estimates, 
  bootciestimates = bootci_estimates
  )

  
  