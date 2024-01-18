#################
# load R package
#################

require(psych)
require(ftsa)
require(tseries)
require(sandwich)
require(Compositional)
require(xtable)

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

#####################
# Interval forecasts
#####################

# object: data matrix
# alpha_val: alpha tuning parameter
# ncomp_tuning: tuning parameter used for selecting the number of components
# fh: forecast horizon
# fore_method: forecast method
# B, K: number of bootstrap samples
# alpha: level of significance

alpha_fun_int <- function(object, alpha_val, ncomp_method, ncomp_tuning = 0.001, fh,
                          fore_method = c("ets", "arima", "rwf"), B = 399, K = 1, alpha)
{
    fore_method = match.arg(fore_method)
    object_alpha = alfa(object, a = alpha_val)$aff
    n_year = nrow(object_alpha)
    n_age = ncol(object_alpha)

    SVD_decomp = svd(object_alpha)
    if(ncomp_method == "eigenvalue")
    {
        ncomp = select_K(tau = ncomp_tuning, SVD_decomp$d^2)
    }
    else if(ncomp_method == "fixed")
    {
        ncomp = 6
    }
    else
    {
        warning("The number of components may be determined by eigenvalue ratio or fixed.")
    }
    basis = as.matrix(SVD_decomp$v[,1:ncomp])
    score = t(object_alpha %*% basis)
    resi = t(object_alpha) - basis %*% score

    # determine in-sample forecast error for principal component scores

    olivia = matrix(NA, ncomp, fh)
    for(ij in 1:ncomp)
    {
        if(fore_method == "ets")
        {
            olivia[ij,] = forecast(ets(score[ij,]), h = fh)$mean
        }
        else if(fore_method == "arima")
        {
            olivia[ij,] = forecast(auto.arima(score[ij,]), h = fh)$mean
        }
        else if(fore_method == "rwf")
        {
            olivia[ij,] = rwf(score[ij,], h = fh, drift = TRUE)$mean
        }
        else
        {
            warning("Please specify a forecasting method listed")
        }
    }
    rownames(olivia) = 1:ncomp
    colnames(olivia) = 1:fh
    forerr = matrix(NA, (n_year - ncomp - fh + 1), ncomp)
    for(i in fh:(n_year - ncomp))
    {
        k = i + (ncomp - fh)
        fore = matrix(NA, 1, ncomp)
        for(j in 1:ncomp)
        {
            if(fore_method == "ets")
            {
                fore[,j] = forecast(ets(score[j,1:k]), h = fh)$mean[fh]
            }
            else if(fore_method == "arima")
            {
                fore[,j] = forecast(auto.arima(score[j,1:k]), h = fh)$mean[fh]
            }
            else if(fore_method == "rwf")
            {
                if(k <= 2)
                {
                    fore[,j] = score[j,k]
                }
                else
                {
                    fore[,j] = rwf(score[j,1:k], h = fh, drift = TRUE)$mean[fh]
                }
            }
        }
        forerr[i-fh+1,] = score[,k+fh] - fore
    }

    # bootstrapping residuals
    q = array(NA, dim = c(n_age, B, K, fh))
    for(j in 1:fh)
    {
        for(i in 1:n_age)
        {
            for(k in 1:K)
            {
                q[i,,k,j] = sample(resi[i,], size = B, replace = TRUE)
            }
        }
    }
    rm(i); rm(j); rm(k)

    # bootstrapping PC score errors

    ny = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
        for(i in 1:ncomp)
        {
            ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE)
        }
    }
    rm(i); rm(j)

    # adding the PC score error to the predicted score

    oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
    fo = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
        for(i in 1:B)
        {
            fo[,i,j] = oli[,i,j] + ny[,i,j]
        }
    }
    rm(i); rm(j)

    # construct bootstrapped samples

    pred = array(NA, dim = c(n_age, B, K, fh))
    for(j in 1:fh)
    {
        for(i in 1:B)
        {
            for(k in 1:K)
            {
                pred[,i,k,j] = basis %*% fo[,i,j] + q[,i,k,j]
            }
        }
    }
    rm(i); rm(j); rm(k)

    pred_resize = array(NA, dim = c(n_age, B * K, fh))
    for(j in 1:fh)
    {
        for(i in 1:B)
        {
            pred_resize[, (((i-1)*K+1):(i*K)), ] = pred[,i,,j]
        }
    }
    rm(i); rm(j)

    # transform back

    d_x_t_star_fore = array(NA, dim = c((n_age + 1), B * K, fh))
    for(iw in 1:fh)
    {
        for(ij in 1:(B * K))
        {
            d_x_t_star_fore[,ij,iw] = as.numeric(alfainv(t(pred_resize[,ij,iw]), a = alpha_val))
        }
    }
    rm(iw); rm(ij)
    return(apply(d_x_t_star_fore, c(1, 3), quantile, c((100 - alpha)/200, 1 - (100 - alpha)/200), na.rm = TRUE))
}

###########
# testing
###########

alpha_int_testing <- function(alpha_para, dat, horizon, method_ncomp, criterion, uni_fore_method, PI_level)
{
    n = nrow(dat)

    # construct prediction intervals

    den_fore = array(NA, dim = c(2, ncol(dat), 11 - horizon))
    for(iw in 1:(11 - horizon))
    {
        den_fore[,,iw] = alpha_fun_int(object = dat[1:(n - 11 + iw),], alpha_val = alpha_para,
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

alpha_CPD_female_testing = ilr_CPD_female_testing = eda_CPD_female_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_female[iwk],
                                                      dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    ilr_CPD_female_testing[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                     horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                     uni_fore_method = "arima", PI_level = 80)

    eda_CPD_female_testing[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue", criterion = "CPD",
                                                    uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female = cbind(alpha_CPD_female_testing, ilr_CPD_female_testing, eda_CPD_female_testing)
xtable(rbind(alpha_CPD_AUS_female, colMeans(alpha_CPD_AUS_female)), digits = 4) #


# interval score

alpha_score_female_testing = ilr_score_female_testing = eda_score_female_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_female[iwk],
                                                        dat = female_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                        criterion = "score", uni_fore_method = "arima", PI_level = 80)

    ilr_score_female_testing[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                       horizon = iwk, method_ncomp = "eigenvalue",
                                                       criterion = "score", uni_fore_method = "arima", PI_level = 80)

    eda_score_female_testing[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                      horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female = cbind(alpha_score_female_testing, ilr_score_female_testing, eda_score_female_testing)
xtable(rbind(alpha_score_AUS_female, colMeans(alpha_score_AUS_female)), digits = 4) # 0.0030 0.0034 0.0075

## male (ARIMA)

# CPD

alpha_CPD_male_testing = ilr_CPD_male_testing = eda_CPD_male_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_male[iwk],
                                                    dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    ilr_CPD_male_testing[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                   horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                   uni_fore_method = "arima", PI_level = 80)

    eda_CPD_male_testing[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                  horizon = iwk, criterion = "CPD", method_ncomp = "eigenvalue",
                                                  uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_male = cbind(alpha_CPD_male_testing, ilr_CPD_male_testing, eda_CPD_male_testing)
xtable(rbind(alpha_CPD_AUS_male, colMeans(alpha_CPD_AUS_male)), digits = 4) # 0.1344 0.1113 0.2552

# interval score

alpha_score_male_testing = ilr_score_male_testing = eda_score_male_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_male[iwk],
                                                      dat = male_pop, horizon = iwk, method_ncomp = "eigenvalue",
                                                      criterion = "score", uni_fore_method = "arima", PI_level = 80)

    ilr_score_male_testing[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                     horizon = iwk, method_ncomp = "eigenvalue",
                                                     criterion = "score", uni_fore_method = "arima", PI_level = 80)

    eda_score_male_testing[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                    horizon = iwk, method_ncomp = "eigenvalue",
                                                    criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male = cbind(alpha_score_male_testing, ilr_score_male_testing, eda_score_male_testing)
xtable(rbind(alpha_score_AUS_male, colMeans(alpha_score_AUS_male)), digits = 4) # 0.0053 0.0052 0.0067

############
# ncomp = 6
############

## female (ARIMA)

# CPD

alpha_CPD_female_testing_ncomp_6 = ilr_CPD_female_testing_ncomp_6 = eda_CPD_female_testing_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_female_ncomp_6[iwk],
                                                            dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    ilr_CPD_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    eda_CPD_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_CPD_AUS_female_ncomp_6 = cbind(alpha_CPD_female_testing_ncomp_6, ilr_CPD_female_testing_ncomp_6, eda_CPD_female_testing_ncomp_6)
xtable(rbind(alpha_CPD_AUS_female_ncomp_6, colMeans(alpha_CPD_AUS_female_ncomp_6)), digits = 4) # 0.0472 0.1062 0.0597

# interval score

alpha_score_female_testing_ncomp_6 = ilr_score_female_testing_ncomp_6 = eda_score_female_testing_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_female_ncomp_6[iwk],
                                                                dat = female_pop, horizon = iwk, method_ncomp = "fixed",
                                                                criterion = "score", uni_fore_method = "arima", PI_level = 80)

    ilr_score_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 0, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)

    eda_score_female_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 1, dat = female_pop,
                                                              horizon = iwk, method_ncomp = "fixed",
                                                              criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_female_ncomp_6 = cbind(alpha_score_female_testing_ncomp_6, ilr_score_female_testing_ncomp_6, eda_score_female_testing_ncomp_6)
xtable(rbind(alpha_score_AUS_female_ncomp_6, colMeans(alpha_score_AUS_female_ncomp_6)), digits = 4) # 0.0023 0.0026 0.0034

## male (ARIMA)

# CPD

alpha_CPD_male_testing_ncomp_6 = ilr_CPD_male_testing_ncomp_6 = eda_CPD_male_testing_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_CPD_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = alpha_para_CPD_AUS_male_ncomp_6[iwk],
                                                            dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    ilr_CPD_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)

    eda_CPD_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "CPD", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)

}

alpha_CPD_AUS_male_ncomp_6 = cbind(alpha_CPD_male_testing_ncomp_6, ilr_CPD_male_testing_ncomp_6, eda_CPD_male_testing_ncomp_6)
xtable(rbind(alpha_CPD_AUS_male_ncomp_6, colMeans(alpha_CPD_AUS_male_ncomp_6)), digits = 4) # 0.1111 0.0724 0.1221

# interval score

alpha_score_male_testing_ncomp_6 = ilr_score_male_testing_ncomp_6 = eda_score_male_testing_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    alpha_score_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = alpha_para_score_AUS_male_ncomp_6[iwk],
                                                            dat = male_pop, horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)

    ilr_score_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 0, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)

    eda_score_male_testing_ncomp_6[iwk] = alpha_int_testing(alpha_para = 1, dat = male_pop,
                                                            horizon = iwk, method_ncomp = "fixed",
                                                            criterion = "score", uni_fore_method = "arima", PI_level = 80)
    print(iwk); rm(iwk)
}

alpha_score_AUS_male_ncomp_6 = cbind(alpha_score_male_testing_ncomp_6, ilr_score_male_testing_ncomp_6, eda_score_male_testing_ncomp_6)
xtable(rbind(alpha_score_AUS_male_ncomp_6, colMeans(alpha_score_AUS_male_ncomp_6)), digits = 4) # 0.0033 0.0030 0.0047
