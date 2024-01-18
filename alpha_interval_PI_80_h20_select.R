########################################################
# Training set 1:(n-41)
# Use validation set (n-40):(n-20) to evaluate accuracy
# Testing set (n-19):n to evaluate accuracy
########################################################

# CPD optim function

alpha_int_cpd_para_select_h20 <- function(alpha_para, dat, horizon, method_ncomp, uni_fore_method, PI_level)
{
    n = nrow(dat)
    
    # construct prediction intervals
    
    den_fore = array(NA, dim = c(2, ncol(dat), 21 - horizon))
    for(iw in 1:(21 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[1:(n - 41 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }
    
    # construct interval scores
    
    int_err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        int_err[ik] = cpd(holdout = dat[(n - 41 + horizon + ik),], lb = den_fore[1,,ik],
                          ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        rm(ik)
    }
    return(mean(int_err))
}

# interval score optim function

alpha_int_score_para_select_h20 <- function(alpha_para, dat, horizon, method_ncomp, uni_fore_method, PI_level)
{
    n = nrow(dat)
    
    # construct prediction intervals
    
    den_fore = array(NA, dim = c(2, ncol(dat), 21 - horizon))
    for(iw in 1:(21 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[1:(n - 41 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }
    
    # construct interval scores
    
    int_err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        int_err[ik] = interval_score(holdout = dat[(n - 41 + horizon + ik),], lb = den_fore[1,,ik],
                                     ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        rm(ik)
    }
    return(mean(int_err))
}

##########
## female 
##########

## eigenvalue ratio criterion

# CPD

alpha_para_CPD_AUS_female_h20 = alpha_para_CPD_AUS_female_value_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_h20, interval = c(0, 1),
                            dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                            uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_female_h20[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_female_h20 = alpha_para_score_AUS_female_value_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_h20, interval = c(0, 1),
                          dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_female_h20[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## fixed (K = 6)

# CPD

alpha_para_CPD_AUS_female_ncomp_6_h20 = alpha_para_CPD_AUS_female_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_h20, interval = c(0, 1),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_female_ncomp_6_h20[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_ncomp_6_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_female_ncomp_6_h20 = alpha_para_score_AUS_female_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_h20, interval = c(0, 1),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_female_ncomp_6_h20[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_ncomp_6_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

#######
# male 
#######

## eigenvalue

# CPD

alpha_para_CPD_AUS_male_h20 = alpha_para_CPD_AUS_male_value_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_h20, interval = c(0, 1),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_male_h20[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_male_h20 = alpha_para_score_AUS_male_value_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_h20, interval = c(0, 1),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_male_h20[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## fixed

# CPD

alpha_para_CPD_AUS_male_ncomp_6_h20 = alpha_para_CPD_AUS_male_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_h20, interval = c(0, 1),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_male_ncomp_6_h20[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_ncomp_6_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_male_ncomp_6_h20 = alpha_para_score_AUS_male_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_h20, interval = c(0, 1),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_male_ncomp_6_h20[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_ncomp_6_h20[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}
