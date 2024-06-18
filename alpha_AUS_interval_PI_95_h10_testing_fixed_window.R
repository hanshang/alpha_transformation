#################
# load R package
#################

require(psych)
require(ftsa)
require(tseries)
require(sandwich)
require(Compositional)

###################################
# Function for computing the score
###################################

# l: lower bound
# u: upper bound
# x: actual holdout data
# alpha: level of significance alpha = 0.2

cpd <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    cover = 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1)))/length(holdout)
    cpd = abs(cover - (1 - alpha))
    rm(lb_ind); rm(ub_ind); rm(cover)
    return(cpd)
}

interval_score <- function(holdout, lb, ub, alpha)
{
    lb_ind = ifelse(holdout < lb, 1, 0)
    ub_ind = ifelse(holdout > ub, 1, 0)
    score = (ub - lb) + 2/alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
    rm(lb_ind); rm(ub_ind)
    return(mean(score))
}

#############
# eigenvalue
#############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_PI_95_fixed_window = ilr_CPD_female_testing_PI_95_fixed_window = eda_CPD_female_testing_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_female_PI_95_fixed_window[iwk],
                                                            dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    ilr_CPD_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                          horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                          uni_fore_method = "arima", PI_level = 95)
    
    eda_CPD_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                          horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                          uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_PI_95_fixed_window = cbind(alpha_CPD_female_testing_PI_95_fixed_window, 
                                                ilr_CPD_female_testing_PI_95_fixed_window, 
                                                eda_CPD_female_testing_PI_95_fixed_window)
colnames(alpha_CPD_AUS_female_PI_95_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_female_testing_PI_95_fixed_window = ilr_score_female_testing_PI_95_fixed_window = eda_score_female_testing_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_female_PI_95_fixed_window[iwk],
                                                                                  dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                                                  criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    ilr_score_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                                                horizon = iwk, method_ncomp = "eigenvalue",
                                                                                criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    eda_score_female_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                                                horizon = iwk, method_ncomp = "eigenvalue",
                                                                                criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_PI_95_fixed_window = cbind(alpha_score_female_testing_PI_95_fixed_window,
                                                  ilr_score_female_testing_PI_95_fixed_window,
                                                  eda_score_female_testing_PI_95_fixed_window)
colnames(alpha_score_AUS_female_PI_95_fixed_window) = c("alpha", "ilr", "eda")

## male (ARIMA)

# CPD

alpha_CPD_male_testing_PI_95_fixed_window = ilr_CPD_male_testing_PI_95_fixed_window = eda_CPD_male_testing_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_male_PI_95_fixed_window[iwk],
                                                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                          criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    ilr_CPD_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                        horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                        uni_fore_method = "arima", PI_level = 95)
    
    eda_CPD_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                        horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                        uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_PI_95_fixed_window = cbind(alpha_CPD_male_testing_PI_95_fixed_window,
                                              ilr_CPD_male_testing_PI_95_fixed_window,
                                              eda_CPD_male_testing_PI_95_fixed_window)
colnames(alpha_CPD_AUS_male_PI_95_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_male_testing_PI_95_fixed_window = ilr_score_male_testing_PI_95_fixed_window = eda_score_male_testing_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_male_PI_95_fixed_window[iwk],
                                                            dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    ilr_score_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                          horizon = iwk, method_ncomp = "eigenvalue",
                                                          criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    eda_score_male_testing_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                          horizon = iwk, method_ncomp = "eigenvalue",
                                                          criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_PI_95_fixed_window = cbind(alpha_score_male_testing_PI_95_fixed_window,
                                                ilr_score_male_testing_PI_95_fixed_window,
                                                eda_score_male_testing_PI_95_fixed_window)
colnames(alpha_score_AUS_male_PI_95_fixed_window) = c("alpha", "ilr", "eda")

############
# ncomp = 6
############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_ncomp_6_PI_95_fixed_window = ilr_CPD_female_testing_ncomp_6_PI_95_fixed_window = eda_CPD_female_testing_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_female_ncomp_6_PI_95_fixed_window[iwk],
                                                                    dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    ilr_CPD_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    eda_CPD_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_ncomp_6_PI_95_fixed_window = cbind(alpha_CPD_female_testing_ncomp_6_PI_95_fixed_window,
                                                        ilr_CPD_female_testing_ncomp_6_PI_95_fixed_window,
                                                        eda_CPD_female_testing_ncomp_6_PI_95_fixed_window)
colnames(alpha_CPD_AUS_female_ncomp_6_PI_95_fixed_window) = c("alpha", "ilr", "eda")
  
# interval score

alpha_score_female_testing_ncomp_6_PI_95_fixed_window = ilr_score_female_testing_ncomp_6_PI_95_fixed_window = eda_score_female_testing_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_female_ncomp_6_PI_95_fixed_window[iwk],
                                                                      dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    ilr_score_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    eda_score_female_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_ncomp_6_PI_95_fixed_window = cbind(alpha_score_female_testing_ncomp_6_PI_95_fixed_window,
                                                         ilr_score_female_testing_ncomp_6_PI_95_fixed_window,
                                                         eda_score_female_testing_ncomp_6_PI_95_fixed_window)
colnames(alpha_score_AUS_female_ncomp_6_PI_95_fixed_window) = c("alpha", "ilr", "eda")

## male (ARIMA)

# CPD

alpha_CPD_male_testing_ncomp_6_PI_95_fixed_window = ilr_CPD_male_testing_ncomp_6_PI_95_fixed_window = eda_CPD_male_testing_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_male_ncomp_6_PI_95_fixed_window[iwk],
                                                                  dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    ilr_CPD_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                                horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    
    eda_CPD_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                                horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_ncomp_6_PI_95_fixed_window = cbind(alpha_CPD_male_testing_ncomp_6_PI_95_fixed_window, 
                                                      ilr_CPD_male_testing_ncomp_6_PI_95_fixed_window, 
                                                      eda_CPD_male_testing_ncomp_6_PI_95_fixed_window)
colnames(alpha_CPD_AUS_male_ncomp_6_PI_95_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_male_testing_ncomp_6_PI_95_fixed_window = ilr_score_male_testing_ncomp_6_PI_95_fixed_window = eda_score_male_testing_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_male_ncomp_6_PI_95_fixed_window[iwk],
                                                                    dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    ilr_score_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "score", uni_fore_method = "arima", PI_level = 95)
    
    eda_score_male_testing_ncomp_6_PI_95_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_ncomp_6_PI_95_fixed_window = cbind(alpha_score_male_testing_ncomp_6_PI_95_fixed_window,
                                                        ilr_score_male_testing_ncomp_6_PI_95_fixed_window,
                                                        eda_score_male_testing_ncomp_6_PI_95_fixed_window)
colnames(alpha_score_AUS_male_ncomp_6_PI_95_fixed_window) = c("alpha", "ilr", "eda")

# summary

xtable(rbind(c(colMeans(alpha_CPD_AUS_female_PI_95_fixed_window), colMeans(alpha_CPD_AUS_female_ncomp_6_PI_95_fixed_window)),
             c(colMeans(alpha_score_AUS_female_PI_95_fixed_window), colMeans(alpha_score_AUS_female_ncomp_6_PI_95_fixed_window)),
            
             c(colMeans(alpha_CPD_AUS_male_PI_95_fixed_window), colMeans(alpha_CPD_AUS_male_ncomp_6_PI_95_fixed_window)),
             c(colMeans(alpha_score_AUS_male_PI_95_fixed_window), colMeans(alpha_score_AUS_male_ncomp_6_PI_95_fixed_window))), digits = 4)
      
