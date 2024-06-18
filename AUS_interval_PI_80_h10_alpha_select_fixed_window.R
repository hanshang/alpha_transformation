########################
# rolling-window scheme
########################

## CPD optim function

alpha_int_cpd_para_select_fixed_window <- function(alpha_para, dat, horizon, method_ncomp, uni_fore_method,
                                                   PI_level)
{
    n = nrow(dat)
    
    # construct prediction intervals
    
    den_fore = array(NA, dim = c(2, ncol(dat), 11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[iw:(n - 21 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }
    
    # construct interval forecast errors
    
    int_err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        int_err[ik] = cpd(holdout = dat[(n - 21 + horizon + ik),], lb = den_fore[1,,ik],
                          ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        rm(ik)
    }
    return(mean(int_err))
}

## interval score optim function

alpha_int_score_para_select_fixed_window <- function(alpha_para, dat, horizon, method_ncomp, uni_fore_method, PI_level)
{
    n = nrow(dat)
    
    # construct prediction intervals
    
    den_fore = array(NA, dim = c(2, ncol(dat), 11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[iw:(n - 21 + iw),], alpha_val = alpha_para,
                                       ncomp_method = method_ncomp, fh = horizon,
                                       fore_method = uni_fore_method, alpha = PI_level)[,,horizon]
        rm(iw)
    }
    
    # construct interval forecast errors
    
    int_err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        int_err[ik] = interval_score(holdout = dat[(n - 21 + horizon + ik),], lb = den_fore[1,,ik],
                                     ub = den_fore[2,,ik], alpha = (100 - PI_level)/100)
        rm(ik)
    }
    return(mean(int_err))
}


##################
## female (ARIMA)
##################

## eigenvalue ratio criterion

# CPD

alpha_para_CPD_AUS_female_fixed_window = alpha_para_CPD_AUS_female_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_female_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_female_fixed_window = alpha_para_score_AUS_female_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_female_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## K = 6

# CPD

alpha_para_CPD_AUS_female_ncomp_6_fixed_window = alpha_para_CPD_AUS_female_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_female_ncomp_6_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_ncomp_6_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_female_ncomp_6_fixed_window = alpha_para_score_AUS_female_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_female_ncomp_6_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_ncomp_6_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}
 
################
## male (ARIMA)
################

## eigenvalue

# CPD

alpha_para_CPD_AUS_male_fixed_window = alpha_para_CPD_AUS_male_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_male_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_male_fixed_window = alpha_para_score_AUS_male_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.45),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_male_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## K = 6

# CPD

alpha_para_CPD_AUS_male_ncomp_6_fixed_window = alpha_para_CPD_AUS_male_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_CPD_AUS_male_ncomp_6_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_ncomp_6_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_male_ncomp_6_fixed_window = alpha_para_score_AUS_male_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 80)
    alpha_para_score_AUS_male_ncomp_6_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_ncomp_6_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

