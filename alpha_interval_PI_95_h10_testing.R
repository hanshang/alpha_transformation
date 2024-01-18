#############
# eigenvalue
#############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_PI_95 = ilr_CPD_female_testing_PI_95 = eda_CPD_female_testing_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_female_PI_95[iwk],
                                                      dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    ilr_CPD_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                     horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                     uni_fore_method = "arima", PI_level = 95)

    eda_CPD_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                    uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_PI_95 = cbind(alpha_CPD_female_testing_PI_95, ilr_CPD_female_testing_PI_95, eda_CPD_female_testing_PI_95)

# interval score

alpha_score_female_testing_PI_95 = ilr_score_female_testing_PI_95 = eda_score_female_testing_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_female_PI_95[iwk],
                                                        dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                        criterion = "score", uni_fore_method = "arima", PI_level = 95)

    ilr_score_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                       horizon = iwk, method_ncomp = "eigenvalue",
                                                       criterion = "score", uni_fore_method = "arima", PI_level = 95)

    eda_score_female_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                      horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_PI_95 = cbind(alpha_score_female_testing_PI_95, ilr_score_female_testing_PI_95, eda_score_female_testing_PI_95)

## male (ARIMA)

# CPD

alpha_CPD_male_testing_PI_95 = ilr_CPD_male_testing_PI_95 = eda_CPD_male_testing_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_male_PI_95[iwk],
                                                    dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    ilr_CPD_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                   horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                   uni_fore_method = "arima", PI_level = 95)

    eda_CPD_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                  horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                  uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_PI_95 = cbind(alpha_CPD_male_testing_PI_95, ilr_CPD_male_testing_PI_95, eda_CPD_male_testing_PI_95)

# interval score

alpha_score_male_testing_PI_95 = ilr_score_male_testing_PI_95 = eda_score_male_testing_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_male_PI_95[iwk],
                                                      dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 95)

    ilr_score_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                     horizon = iwk, method_ncomp = "eigenvalue",
                                                     criterion = "score", uni_fore_method = "arima", PI_level = 95)

    eda_score_male_testing_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_PI_95 = cbind(alpha_score_male_testing_PI_95, ilr_score_male_testing_PI_95, eda_score_male_testing_PI_95)

############
# ncomp = 6
############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_ncomp_6_PI_95 = ilr_CPD_female_testing_ncomp_6_PI_95 = eda_CPD_female_testing_ncomp_6_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_female_ncomp_6_PI_95[iwk],
                                                            dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    ilr_CPD_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    eda_CPD_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_ncomp_6_PI_95 = cbind(alpha_CPD_female_testing_ncomp_6_PI_95, ilr_CPD_female_testing_ncomp_6_PI_95, eda_CPD_female_testing_ncomp_6_PI_95)

# interval score

alpha_score_female_testing_ncomp_6_PI_95 = ilr_score_female_testing_ncomp_6_PI_95 = eda_score_female_testing_ncomp_6_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_female_ncomp_6_PI_95[iwk],
                                                                      dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                      criterion = "score", uni_fore_method = "arima", PI_level = 95)

    ilr_score_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)

    eda_score_female_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                                    horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_ncomp_6_PI_95 = cbind(alpha_score_female_testing_ncomp_6_PI_95, ilr_score_female_testing_ncomp_6_PI_95, eda_score_female_testing_ncomp_6_PI_95)

## male (ARIMA)

# CPD

alpha_CPD_male_testing_ncomp_6_PI_95 = ilr_CPD_male_testing_ncomp_6_PI_95 = eda_CPD_male_testing_ncomp_6_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_male_ncomp_6_PI_95[iwk],
                                                                  dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    ilr_CPD_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                                horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "CPD", uni_fore_method = "arima", PI_level = 95)

    eda_CPD_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                                horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "CPD", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male_ncomp_6_PI_95 = cbind(alpha_CPD_male_testing_ncomp_6_PI_95, ilr_CPD_male_testing_ncomp_6_PI_95, eda_CPD_male_testing_ncomp_6_PI_95)

# interval score

alpha_score_male_testing_ncomp_6_PI_95 = ilr_score_male_testing_ncomp_6_PI_95 = eda_score_male_testing_ncomp_6_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_male_ncomp_6_PI_95[iwk],
                                                                    dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                                    criterion = "score", uni_fore_method = "arima", PI_level = 95)

    ilr_score_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "score", uni_fore_method = "arima", PI_level = 95)

    eda_score_male_testing_ncomp_6_PI_95[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                                  horizon = iwk, method_ncomp = "fixed",
                                                                  criterion = "score", uni_fore_method = "arima", PI_level = 95)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_ncomp_6_PI_95 = cbind(alpha_score_male_testing_ncomp_6_PI_95, ilr_score_male_testing_ncomp_6_PI_95, eda_score_male_testing_ncomp_6_PI_95)
