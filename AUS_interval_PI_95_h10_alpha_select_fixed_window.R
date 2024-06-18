############
# PI = 95%
############

## eigenvalue

# CPD

alpha_para_CPD_AUS_female_PI_95_fixed_window = alpha_para_CPD_AUS_female_value_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_CPD_AUS_female_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_female_PI_95_fixed_window = alpha_para_score_AUS_female_value_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_score_AUS_female_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## K = 6

alpha_para_CPD_AUS_female_ncomp_6_PI_95_fixed_window = alpha_para_CPD_AUS_female_value_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_CPD_AUS_female_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_female_value_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

alpha_para_score_AUS_female_ncomp_6_PI_95_fixed_window = alpha_para_score_AUS_female_value_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_score_AUS_female_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_female_value_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

### male

## eigenvalue

# CPD

alpha_para_CPD_AUS_male_PI_95_fixed_window = alpha_para_CPD_AUS_male_value_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_CPD_AUS_male_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

# interval score

alpha_para_score_AUS_male_PI_95_fixed_window = alpha_para_score_AUS_male_value_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_score_AUS_male_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

## K = 6

alpha_para_CPD_AUS_male_ncomp_6_PI_95_fixed_window = alpha_para_CPD_AUS_male_value_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_cpd_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_CPD_AUS_male_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_CPD_AUS_male_value_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}

alpha_para_score_AUS_male_ncomp_6_PI_95_fixed_window = alpha_para_score_AUS_male_value_ncomp_6_PI_95_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    optim_fun <- optimise(f = alpha_int_score_para_select_fixed_window, interval = c(0, 0.5),
                          dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                          uni_fore_method = "arima", PI_level = 95)
    alpha_para_score_AUS_male_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$minimum
    alpha_para_score_AUS_male_value_ncomp_6_PI_95_fixed_window[iwk] = optim_fun$objective
    print(iwk); rm(iwk); rm(optim_fun)
}



