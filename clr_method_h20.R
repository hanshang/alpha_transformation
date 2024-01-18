########################
# CoDa testing (h = 20)
########################

CoDa_testing_h20 <- function(dat, horizon, criterion, selecting_ncomp, uni_method)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 21 - horizon)
    for(iw in 1:(21 - horizon))
    {
        den_fore[,iw] = R_square_fit(dat = dat[1:(n - 21 + iw),], fh = horizon, ncomp_selection = selecting_ncomp,
                                    forecasting_method = uni_method)$fore_count[,horizon]
    }

    err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        data_c = cbind(True = dat[(n - 21 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
        data_c = data_c[!is.na(data_c[,1]),]
        if(any(which(!is.finite(data_c))))
        {
            err[ik] = "NA"
        }
        else
        {
            if(criterion == "KL_div")
            {
                err[ik] = mean(as.numeric(KLdiv(data_c, eps = 1e-16))[2:3])
            }
            else if(criterion == "JS_div_simple")
            {
                M = rowMeans(data_c)
                P_M = cbind(data_c[,1], M)
                E_M = cbind(as.numeric(data_c[,2]), M)
                colnames(E_M) = colnames(P_M) = c("True", "M")
                err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
                rm(M); rm(P_M); rm(E_M)
            }
            else if(criterion == "JS_div_geo")
            {
                M = apply(data_c, 1, geometric.mean)
                P_M = cbind(data_c[,1], M)
                E_M = cbind(as.numeric(data_c[,2]), M)
                colnames(E_M) = colnames(P_M) = c("True", "M")
                err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
                rm(M); rm(P_M); rm(E_M)
            }
            else
            {
                warning("Please specify a criterion from the list.")
            }
        }
    }
    rm(n); rm(den_fore)
    return(mean(as.numeric(err), na.rm = TRUE))
}

###################
# eigenvalue ratio
###################

CoDa_KL_div_female_testing_h20 = CoDa_KL_div_male_testing_h20 =
CoDa_JS_div_simple_female_testing_h20 = CoDa_JS_div_simple_male_testing_h20 =
CoDa_JS_div_geo_female_testing_h20 = CoDa_JS_div_geo_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    # KL_div

    CoDa_KL_div_female_testing_h20[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                   selecting_ncomp = "eigen_value_ratio",
                                                   uni_method = "ARIMA")

    CoDa_KL_div_male_testing_h20[iwk]   = CoDa_testing_h20(dat = male_pop,   horizon = iwk, criterion = "KL_div",
                                                   selecting_ncomp = "eigen_value_ratio",
                                                   uni_method = "ARIMA")

    # JS_div_simple

    CoDa_JS_div_simple_female_testing_h20[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "eigen_value_ratio",
                                                          uni_method = "ARIMA")

    CoDa_JS_div_simple_male_testing_h20[iwk]   = CoDa_testing_h20(dat = male_pop,   horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "eigen_value_ratio",
                                                          uni_method = "ARIMA")

    # JS_div_geo

    CoDa_JS_div_geo_female_testing_h20[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "eigen_value_ratio",
                                                       uni_method = "ARIMA")

    CoDa_JS_div_geo_male_testing_h20[iwk]   = CoDa_testing_h20(dat = male_pop,   horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "eigen_value_ratio",
                                                       uni_method = "ARIMA")
    print(iwk); rm(iwk)
}

CoDa_results_testing_h20 = cbind(CoDa_KL_div_female_testing_h20, CoDa_JS_div_simple_female_testing_h20, CoDa_JS_div_geo_female_testing_h20,
                                 CoDa_KL_div_male_testing_h20,   CoDa_JS_div_simple_male_testing_h20,   CoDa_JS_div_geo_male_testing_h20)
colnames(CoDa_results_testing_h20) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                       "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(CoDa_results_testing_h20) = 1:20

##################
# CoDa_LC (K = 6)
##################

CoDa_KL_div_female_testing_h20_ncomp_6 = CoDa_KL_div_male_testing_h20_ncomp_6 =
CoDa_JS_div_simple_female_testing_h20_ncomp_6 = CoDa_JS_div_simple_male_testing_h20_ncomp_6 =
CoDa_JS_div_geo_female_testing_h20_ncomp_6 = CoDa_JS_div_geo_male_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    # KL_div

    CoDa_KL_div_female_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                                   selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_KL_div_male_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = male_pop, horizon = iwk, criterion = "KL_div",
                                                                 selecting_ncomp = "fixed", uni_method = "ARIMA")

    # JS_div_simple

    CoDa_JS_div_simple_female_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "JS_div_simple",
                                                                          selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_JS_div_simple_male_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = male_pop, horizon = iwk, criterion = "JS_div_simple",
                                                                        selecting_ncomp = "fixed", uni_method = "ARIMA")

    # JS_div_geo

    CoDa_JS_div_geo_female_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = female_pop, horizon = iwk, criterion = "JS_div_geo",
                                                                       selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_JS_div_geo_male_testing_h20_ncomp_6[iwk] = CoDa_testing_h20(dat = male_pop, horizon = iwk, criterion = "JS_div_geo",
                                                                     selecting_ncomp = "fixed", uni_method = "ARIMA")
    print(iwk); rm(iwk)
}

CoDa_results_testing_h20_ncomp_6 = cbind(CoDa_KL_div_female_testing_h20_ncomp_6, CoDa_JS_div_simple_female_testing_h20_ncomp_6, CoDa_JS_div_geo_female_testing_h20_ncomp_6,
                                         CoDa_KL_div_male_testing_h20_ncomp_6,   CoDa_JS_div_simple_male_testing_h20_ncomp_6,   CoDa_JS_div_geo_male_testing_h20_ncomp_6)
colnames(CoDa_results_testing_h20_ncomp_6) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                               "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(CoDa_results_testing_h20_ncomp_6) = 1:20
