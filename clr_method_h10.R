#############
# CoDa model
#############

# dat: (n by p) data matrix
# fh: forecast horizon
# ncomp_tuning: tuning parameter for selecting the number of components
# forecasting_method: univariate time-series method

R_square_fit <- function(dat, fh, ncomp_selection, ncomp_tuning = 0.001,
                         forecasting_method)
{
    n_year = nrow(dat)
    n_age = ncol(dat)

    # standardize life table death to sum to 1

    dat_center = sweep(dat, 1, apply(dat, 1, sum), "/")

    alpha_x = vector("numeric", n_age)
    for(ik in 1:n_age)
    {
        alpha_x[ik] = geometric.mean(dat_center[,ik])
    }

    f_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
        f_x_t[ik,] = (dat[ik,]/alpha_x)/sum(dat[ik,]/alpha_x)
    }

    g_t = vector("numeric", n_year)
    h_x_t = matrix(NA, n_year, n_age)
    for(ik in 1:n_year)
    {
        g_t[ik] = geometric.mean(f_x_t[ik,])
        h_x_t[ik,] = log(f_x_t[ik,]/g_t[ik])
    }

    SVD_decomp = svd(h_x_t) # 1st component only
    if(ncomp_selection == "one")
    {
        ncomp = 1
    }
    else if(ncomp_selection == "eigen_value_ratio")
    {
        ncomp = select_K(tau = ncomp_tuning, SVD_decomp$d^2)
    }
    else if(ncomp_selection == "fixed")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of components should be correctly selected.")
    }
    basis = SVD_decomp$v[,1:ncomp]
    score = t(basis) %*% t(h_x_t)
    recon = basis %*% score

    # reconstruction (model in-sample fitting)

    f_x_t_star_recon = d_x_t_star_recon = matrix(NA, n_age, n_year)
    for(ik in 1:n_year)
    {
        f_x_t_star_recon[,ik] = exp(recon[,ik])/sum(exp(recon[,ik]))
        d_x_t_star_recon[,ik] = (f_x_t_star_recon[,ik] * alpha_x)/sum(f_x_t_star_recon[,ik] * alpha_x)
    }
    R2 = 1 - sum((t(d_x_t_star_recon) * 100000 - dat)^2)/sum((dat - colMeans(dat))^2)

    # forecasts of principal component scores

    score_fore = matrix(NA, ncomp, fh)
    for(ik in 1:ncomp)
    {
        if(forecasting_method == "RWF_no_drift")
        {
            score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = FALSE)$mean
        }
        else if(forecasting_method == "RWF_drift")
        {
            score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = TRUE)$mean
        }
        else if(forecasting_method == "ETS")
        {
            score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean
        }
        else if(forecasting_method == "ARIMA")
        {
            score_fore[ik,] = forecast(auto.arima(as.numeric(score[ik,])), h = fh)$mean
        }
        else
        {
            warning("Univariate time series method is not among the list.")
        }
    }

    # obtain forecasts in real-valued space

    fore_val = basis %*% score_fore

    # back-transformation

    f_x_t_star_fore = d_x_t_star_fore = matrix(NA, n_age, fh)
    for(ik in 1:fh)
    {
        f_x_t_star_fore[,ik] = exp(fore_val[,ik])/sum(exp(fore_val[,ik]))
        d_x_t_star_fore[,ik] = (f_x_t_star_fore[,ik] * alpha_x)/sum((f_x_t_star_fore[,ik] * alpha_x))
    }
    return(list(R2 = R2, recon = d_x_t_star_recon, fore_count = d_x_t_star_fore,
                alpha_x = alpha_x, h_x_t = h_x_t, basis_fore = basis, score_fore = score_fore,
                score = score, ncomp = ncomp))
}

########################
# CoDa testing (h = 10)
########################

CoDa_testing <- function(dat, horizon, criterion, selecting_ncomp, uni_method)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 11 - horizon)
    for(iw in 1:(11 - horizon))
    {
        den_fore[,iw] = R_square_fit(dat = dat[1:(n - 11 + iw),], fh = horizon, ncomp_selection = selecting_ncomp,
                                     forecasting_method = uni_method)$fore_count[,horizon]
    }

    err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        data_c = cbind(True = dat[(n - 11 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
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

CoDa_KL_div_female_testing = CoDa_KL_div_male_testing =
CoDa_JS_div_simple_female_testing = CoDa_JS_div_simple_male_testing =
CoDa_JS_div_geo_female_testing = CoDa_JS_div_geo_male_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    # KL_div

    CoDa_KL_div_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                   selecting_ncomp = "eigen_value_ratio",
                                                   uni_method = "ARIMA")

    CoDa_KL_div_male_testing[iwk]   = CoDa_testing(dat = male_pop,   horizon = iwk, criterion = "KL_div",
                                                   selecting_ncomp = "eigen_value_ratio",
                                                   uni_method = "ARIMA")

    # JS_div_simple

    CoDa_JS_div_simple_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "eigen_value_ratio",
                                                          uni_method = "ARIMA")

    CoDa_JS_div_simple_male_testing[iwk]   = CoDa_testing(dat = male_pop,   horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "eigen_value_ratio",
                                                          uni_method = "ARIMA")

    # JS_div_geo

    CoDa_JS_div_geo_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "eigen_value_ratio",
                                                       uni_method = "ARIMA")

    CoDa_JS_div_geo_male_testing[iwk]   = CoDa_testing(dat = male_pop,   horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "eigen_value_ratio",
                                                       uni_method = "ARIMA")
    print(iwk); rm(iwk)
}

CoDa_results_testing = cbind(CoDa_KL_div_female_testing, CoDa_JS_div_simple_female_testing, CoDa_JS_div_geo_female_testing,
                             CoDa_KL_div_male_testing,   CoDa_JS_div_simple_male_testing,   CoDa_JS_div_geo_male_testing)
colnames(CoDa_results_testing) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                   "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(CoDa_results_testing) = 1:10

##################
# CoDa_LC (K = 1)
##################

CoDa_LC_KL_div_female_testing = CoDa_LC_KL_div_male_testing =
CoDa_LC_JS_div_simple_female_testing = CoDa_LC_JS_div_simple_male_testing =
CoDa_LC_JS_div_geo_female_testing = CoDa_LC_JS_div_geo_male_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    # KL_div

    CoDa_LC_KL_div_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                      selecting_ncomp = "one", uni_method = "RWF_drift")

    CoDa_LC_KL_div_male_testing[iwk]   = CoDa_testing(dat = male_pop, horizon = iwk, criterion = "KL_div",
                                                      selecting_ncomp = "one", uni_method = "RWF_drift")

    # JS_div_simple

    CoDa_JS_div_simple_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "one", uni_method = "RWF_drift")

    CoDa_JS_div_simple_male_testing[iwk]   = CoDa_testing(dat = male_pop,   horizon = iwk, criterion = "JS_div_simple",
                                                          selecting_ncomp = "one", uni_method = "RWF_drift")

    # JS_div_geo

    CoDa_JS_div_geo_female_testing[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "one", uni_method = "RWF_drift")

    CoDa_JS_div_geo_male_testing[iwk]   = CoDa_testing(dat = male_pop,   horizon = iwk, criterion = "JS_div_geo",
                                                       selecting_ncomp = "one", uni_method = "RWF_drift")
    print(iwk); rm(iwk)
}

CoDa_LC_results_testing = cbind(CoDa_LC_KL_div_female_testing, CoDa_JS_div_simple_female_testing, CoDa_JS_div_geo_female_testing,
                                CoDa_LC_KL_div_male_testing,   CoDa_JS_div_simple_male_testing,   CoDa_JS_div_geo_male_testing)
colnames(CoDa_LC_results_testing) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                      "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(CoDa_LC_results_testing) = 1:10

##################
# CoDa_LC (K = 6)
##################

CoDa_KL_div_female_testing_ncomp_6 = CoDa_KL_div_male_testing_ncomp_6 =
CoDa_JS_div_simple_female_testing_ncomp_6 = CoDa_JS_div_simple_male_testing_ncomp_6 =
CoDa_JS_div_geo_female_testing_ncomp_6 = CoDa_JS_div_geo_male_testing_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    # KL_div

    CoDa_KL_div_female_testing_ncomp_6[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                           selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_KL_div_male_testing_ncomp_6[iwk] = CoDa_testing(dat = male_pop, horizon = iwk, criterion = "KL_div",
                                                         selecting_ncomp = "fixed", uni_method = "ARIMA")

    # JS_div_simple

    CoDa_JS_div_simple_female_testing_ncomp_6[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_simple",
                                                                  selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_JS_div_simple_male_testing_ncomp_6[iwk] = CoDa_testing(dat = male_pop, horizon = iwk, criterion = "JS_div_simple",
                                                                selecting_ncomp = "fixed", uni_method = "ARIMA")

    # JS_div_geo

    CoDa_JS_div_geo_female_testing_ncomp_6[iwk] = CoDa_testing(dat = female_pop, horizon = iwk, criterion = "JS_div_geo",
                                                               selecting_ncomp = "fixed", uni_method = "ARIMA")

    CoDa_JS_div_geo_male_testing_ncomp_6[iwk] = CoDa_testing(dat = male_pop, horizon = iwk, criterion = "JS_div_geo",
                                                             selecting_ncomp = "fixed", uni_method = "ARIMA")
    print(iwk); rm(iwk)
}

CoDa_results_testing_ncomp_6 = cbind(CoDa_KL_div_female_testing_ncomp_6, CoDa_JS_div_simple_female_testing_ncomp_6, CoDa_JS_div_geo_female_testing_ncomp_6,
                                     CoDa_KL_div_male_testing_ncomp_6,   CoDa_JS_div_simple_male_testing_ncomp_6,   CoDa_JS_div_geo_male_testing_ncomp_6)
colnames(CoDa_results_testing_ncomp_6) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                           "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(CoDa_results_testing_ncomp_6) = 1:10
