###########
# testing
###########

alpha_int_testing_fixed_window <- function(alpha_para, dat, horizon, method_ncomp, criterion, uni_fore_method, PI_level)
{
    n = nrow(dat)

    # construct prediction intervals

    den_fore = array(NA, dim = c(2, ncol(dat), 11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[iw:(n - 11 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }

    # construct interval forecast errors

    int_err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        if(criterion == "CPD")
        {
            int_err[ik] = cpd(holdout = dat[(n - 11 + horizon + ik),], lb = den_fore[1,,ik],
                                      ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        }
        else if(criterion == "score")
        {
            int_err[ik] = interval_score(holdout = dat[(n - 11 + horizon + ik),], lb = den_fore[1,,ik],
                                         ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        }
        else
        {
            warning("Criterion must be either CPD or interval score")
        }
        rm(ik)
    }
    result_ave = mean(int_err)
    rm(int_err); rm(den_fore); rm(n)
    return(result_ave)
}

#############
# eigenvalue
#############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_fixed_window = ilr_CPD_female_testing_fixed_window = eda_CPD_female_testing_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_female_fixed_window[iwk],
                                                      dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    ilr_CPD_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                     horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                     uni_fore_method = "arima", PI_level = 80)

    eda_CPD_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                    uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_fixed_window = cbind(alpha_CPD_female_testing_fixed_window, ilr_CPD_female_testing_fixed_window, 
                                          eda_CPD_female_testing_fixed_window)
colnames(alpha_CPD_AUS_female_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_female_testing_fixed_window = ilr_score_female_testing_fixed_window = eda_score_female_testing_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_female_fixed_window[iwk],
                                                        dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                        criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                      horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_female_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                      horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_fixed_window = cbind(alpha_score_female_testing_fixed_window, ilr_score_female_testing_fixed_window, 
                                            eda_score_female_testing_fixed_window)
colnames(alpha_score_AUS_female_fixed_window) = c("alpha", "ilr", "eda")

## male (ARIMA)

# CPD

alpha_CPD_male_testing_fixed_window = ilr_CPD_male_testing_fixed_window = eda_CPD_male_testing_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_male_fixed_window[iwk],
                                                    dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                  horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                  uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                  horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                  uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_fixed_window = cbind(alpha_CPD_male_testing_fixed_window, ilr_CPD_male_testing_fixed_window, 
                                        eda_CPD_male_testing_fixed_window)
colnames(alpha_CPD_AUS_male_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_male_testing_fixed_window = ilr_score_male_testing_fixed_window = eda_score_male_testing_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_male_fixed_window[iwk],
                                                      dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_male_testing_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_fixed_window = cbind(alpha_score_male_testing_fixed_window, ilr_score_male_testing_fixed_window, 
                                          eda_score_male_testing_fixed_window)
colnames(alpha_score_AUS_male_fixed_window) = c("alpha", "ilr", "eda")

############
# ncomp = 6
############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_ncomp_6_fixed_window = ilr_CPD_female_testing_ncomp_6_fixed_window = eda_CPD_female_testing_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_female_ncomp_6_fixed_window[iwk],
                                                              dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_ncomp_6_fixed_window = cbind(alpha_CPD_female_testing_ncomp_6_fixed_window,
                                                  ilr_CPD_female_testing_ncomp_6_fixed_window,
                                                  eda_CPD_female_testing_ncomp_6_fixed_window)
colnames(alpha_CPD_AUS_female_ncomp_6_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_female_testing_ncomp_6_fixed_window = ilr_score_female_testing_ncomp_6_fixed_window = eda_score_female_testing_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_female_ncomp_6_fixed_window[iwk],
                                                                dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_female_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_ncomp_6_fixed_window = cbind(alpha_score_female_testing_ncomp_6_fixed_window, 
                                                    ilr_score_female_testing_ncomp_6_fixed_window, 
                                                    eda_score_female_testing_ncomp_6_fixed_window)
colnames(alpha_score_AUS_female_ncomp_6_fixed_window) = c("alpha", "ilr", "eda")

## male (ARIMA)

# CPD

alpha_CPD_male_testing_ncomp_6_fixed_window = ilr_CPD_male_testing_ncomp_6_fixed_window = eda_CPD_male_testing_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_CPD_AUS_male_ncomp_6_fixed_window[iwk],
                                                            dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                          horizon = iwk, method_ncomp = "fixed",
                                                          criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                          horizon = iwk, method_ncomp = "fixed",
                                                          criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_ncomp_6_fixed_window = cbind(alpha_CPD_male_testing_ncomp_6_fixed_window, 
                                                ilr_CPD_male_testing_ncomp_6_fixed_window, 
                                                eda_CPD_male_testing_ncomp_6_fixed_window)
colnames(alpha_CPD_AUS_male_ncomp_6_fixed_window) = c("alpha", "ilr", "eda")

# interval score

alpha_score_male_testing_ncomp_6_fixed_window = ilr_score_male_testing_ncomp_6_fixed_window = 
eda_score_male_testing_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = alpha_para_score_AUS_male_ncomp_6_fixed_window[iwk],
                                                              dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 0, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_male_testing_ncomp_6_fixed_window[iwk] = alpha_int_testing_fixed_window(alpha_para = 1, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_ncomp_6_fixed_window = cbind(alpha_score_male_testing_ncomp_6_fixed_window, 
                                                  ilr_score_male_testing_ncomp_6_fixed_window, 
                                                  eda_score_male_testing_ncomp_6_fixed_window)
colnames(alpha_score_AUS_male_ncomp_6_fixed_window) = c("alpha", "ilr", "eda")

# summary

require(xtable)
xtable(rbind(c(colMeans(alpha_CPD_AUS_female_fixed_window), colMeans(alpha_CPD_AUS_female_ncomp_6_fixed_window)),
            c(colMeans(alpha_score_AUS_female_fixed_window), colMeans(alpha_score_AUS_female_ncomp_6_fixed_window)),
            
            c(colMeans(alpha_CPD_AUS_male_fixed_window), colMeans(alpha_CPD_AUS_male_ncomp_6_fixed_window)),
            c(colMeans(alpha_score_AUS_male_fixed_window), colMeans(alpha_score_AUS_male_ncomp_6_fixed_window))), digits = 4)

