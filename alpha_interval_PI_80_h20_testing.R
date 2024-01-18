###########
# testing
###########

alpha_int_testing_h20 <- function(alpha_para, dat, horizon, method_ncomp, criterion, uni_fore_method, PI_level)
{
    n = nrow(dat)
    
    # construct prediction intervals
    
    den_fore = array(NA, dim = c(2, ncol(dat), 21 - horizon))
    for(iw in 1:(21 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[1:(n - 21 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }
    
    # construct interval forecast errors
    
    int_err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        if(criterion == "CPD")
        {
            int_err[ik] = cpd(holdout = dat[(n - 21 + horizon + ik),], lb = den_fore[1,,ik],
                              ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        }
        else if(criterion == "score")
        {
            int_err[ik] = interval_score(holdout = dat[(n - 21 + horizon + ik),], lb = den_fore[1,,ik],
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

alpha_CPD_female_testing_h20 = ilr_CPD_female_testing_h20 = eda_CPD_female_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_CPD_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_CPD_AUS_female_h20[iwk],
                                                              dat = female_pop, horizon = iwk, 
                                                              method_ncomp = "eigenvalue", criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
            
    ilr_CPD_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                            uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                            uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_h20 = round(colMeans(cbind(alpha_CPD_female_testing_h20, 
                                                ilr_CPD_female_testing_h20, 
                                                eda_CPD_female_testing_h20)), 4)

# interval score

alpha_score_female_testing_h20 = ilr_score_female_testing_h20 = eda_score_female_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_score_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_score_AUS_female_h20[iwk],
                                                                dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                                criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "eigenvalue",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_female_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "eigenvalue",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_h20 = round(colMeans(cbind(alpha_score_female_testing_h20, 
                                                  ilr_score_female_testing_h20, 
                                                  eda_score_female_testing_h20)), 4)

## male (ARIMA)

# CPD

alpha_CPD_male_testing_h20 = ilr_CPD_male_testing_h20 = eda_CPD_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_CPD_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_CPD_AUS_male_h20[iwk],
                                                            dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = male_pop,
                                                          horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                          uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = male_pop,
                                                          horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                          uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_h20 = round(colMeans(cbind(alpha_CPD_male_testing_h20, 
                                              ilr_CPD_male_testing_h20, 
                                              eda_CPD_male_testing_h20)), 4)

# interval score

alpha_score_male_testing_h20 = ilr_score_male_testing_h20 = eda_score_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_score_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_score_AUS_male_h20[iwk],
                                                              dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "eigenvalue",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_male_testing_h20[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "eigenvalue",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_h20 = round(colMeans(cbind(alpha_score_male_testing_h20, 
                                                ilr_score_male_testing_h20, 
                                                eda_score_male_testing_h20)), 4)

########
# K = 6
########

## female (ARIMA)

# CPD

alpha_CPD_female_testing_h20_ncomp_6 = ilr_CPD_female_testing_h20_ncomp_6 = eda_CPD_female_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_CPD_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_CPD_AUS_female_ncomp_6_h20[iwk],
                                                                      dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_h20_ncomp_6 = round(colMeans(cbind(alpha_CPD_female_testing_h20_ncomp_6, 
                                                        ilr_CPD_female_testing_h20_ncomp_6, 
                                                        eda_CPD_female_testing_h20_ncomp_6)), 4)

# interval score

alpha_score_female_testing_h20_ncomp_6 = ilr_score_female_testing_h20_ncomp_6 = eda_score_female_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_score_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_score_AUS_female_ncomp_6_h20[iwk],
                                                                        dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                        criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    ilr_score_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = female_pop,
                                                                      horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_female_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = female_pop,
                                                                      horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_h20_ncomp_6 = round(colMeans(cbind(alpha_score_female_testing_h20_ncomp_6, 
                                                          ilr_score_female_testing_h20_ncomp_6, 
                                                          eda_score_female_testing_h20_ncomp_6)), 4)

## male (ARIMA)

# CPD

alpha_CPD_male_testing_h20_ncomp_6 = ilr_CPD_male_testing_h20_ncomp_6 = eda_CPD_male_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_CPD_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_CPD_AUS_male_ncomp_6_h20[iwk],
                                                                    dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    ilr_CPD_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    
    eda_CPD_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_h20_ncomp_6 = round(colMeans(cbind(alpha_CPD_male_testing_h20_ncomp_6, 
                                                      ilr_CPD_male_testing_h20_ncomp_6, 
                                                      eda_CPD_male_testing_h20_ncomp_6)), 4)

# interval score

alpha_score_male_testing_h20_ncomp_6 = ilr_score_male_testing_h20_ncomp_6 = eda_score_male_testing_h20_ncomp_6 = vector("numeric", 10)
for(iwk in 1:20)
{
    alpha_score_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = alpha_para_score_AUS_male_ncomp_6_h20[iwk],
                                                                      dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
            
    ilr_score_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 0, dat = male_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 80)
    
    eda_score_male_testing_h20_ncomp_6[iwk] = alpha_int_testing_h20(alpha_para = 1, dat = male_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_h20_ncomp_6 = round(colMeans(cbind(alpha_score_male_testing_h20_ncomp_6, 
                                                        ilr_score_male_testing_h20_ncomp_6, 
                                                        eda_score_male_testing_h20_ncomp_6)), 4)
