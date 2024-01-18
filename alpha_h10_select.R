#######################
# alpha transformation
#######################

source("choice_K.R")

# object: life-table death counts
# alpha_val: alpha between 0 and 1
# ncomp_tuning: tuning parameter in the selection of components
# fh: forecast horizon
# Helmert_mat: TRUE or FALSE

alpha_fun <- function(object, alpha_val, ncomp_method, ncomp_tuning = 0.001, fh, Helmert_mat)
{
    object_alpha = alfa(object, a = alpha_val, h = Helmert_mat)$aff
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
        warning("The number of components can be either determined by eigenvalue or fixed.")
    }
    
    basis = as.matrix(SVD_decomp$v[,1:ncomp])
    score = object_alpha %*% basis
    recon = t(basis %*% t(score))
    
    # reconstruction
    
    object_recon = alfainv(recon, a = alpha_val, h = Helmert_mat)
    R2 = round(1 - sum((object_recon - object)^2)/sum((object - colMeans(object))^2), 4)
    
    # forecasts of principal component scores
    
    score_fore = matrix(NA, ncomp, fh)
    for(ik in 1:ncomp)
    {
        score_fore[ik,] = forecast(auto.arima(as.numeric(score[,ik])), h = fh)$mean
    }
    
    # obtain forecasts in real-valued space
    
    fore_val = t(basis %*% score_fore)
    object_fore = alfainv(fore_val, a = alpha_val, h = Helmert_mat)
    return(t(object_fore))
}


########################################################
# Training set 1:(n-21)
# Use validation set (n-20):(n-10) to evaluate accuracy
# Testing set (n-9):n to evaluate accuracy
########################################################

alpha_para_select <- function(alpha_para, dat, horizon, criterion, mat_Helmert, method_ncomp)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 11 - horizon)
    for(iw in 1:(11 - horizon))
    {
        den_fore[,iw] = alpha_fun(object = dat[1:(n - 21 + iw),], alpha_val = alpha_para, 
                                  ncomp_method = method_ncomp, 
                                  fh = horizon, Helmert_mat = mat_Helmert)[,horizon]
    }
    
    # true_dat = replace(dat, which(dat == 0), NA)
    err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
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

##########
## KL_div
##########

## Helmert = TRUE (sum to 0)

# AUS_female

alpha_para_KL_div_AUS_female = alpha_para_KL_div_AUS_female_value = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk,  criterion = "KL_div",
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_female[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_KL_div_AUS_female_ncomp_6 = alpha_para_KL_div_AUS_female_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk, criterion = "KL_div", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_female_ncomp_6[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# AUS_male

alpha_para_KL_div_AUS_male = alpha_para_KL_div_AUS_male_value = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk,  
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_male[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_KL_div_AUS_male_ncomp_6 = alpha_para_KL_div_AUS_male_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk, criterion = "KL_div", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_male_ncomp_6[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

#################
## JS_div_simple
#################

## Helmert = TRUE (sum to 0)

# AUS_female

alpha_para_JS_div_simple_AUS_female = alpha_para_JS_div_simple_AUS_female_value = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_simple_AUS_female[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_simple_AUS_female_ncomp_6 = alpha_para_JS_div_simple_AUS_female_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_simple_AUS_female_ncomp_6[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# AUS_male

alpha_para_JS_div_simple_AUS_male = alpha_para_JS_div_simple_AUS_male_value = vector("numeric", 10)
for(iwk in 8:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk,  criterion = "JS_div_simple",
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_simple_AUS_male[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_simple_AUS_male_ncomp_6 = alpha_para_JS_div_simple_AUS_male_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_simple_AUS_male_ncomp_6[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

#############
# JS_div_geo
#############

## Helmert = TRUE 

# AUS_female

alpha_para_JS_div_geo_AUS_female = alpha_para_JS_div_geo_AUS_female_value = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk,  criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_female[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_geo_AUS_female_ncomp_6 = alpha_para_JS_div_geo_AUS_female_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = female_pop, horizon = iwk, criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_female_ncomp_6[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

## AUS_male

alpha_para_JS_div_geo_AUS_male = alpha_para_JS_div_geo_AUS_male_value = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk,  criterion = "JS_div_geo",
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_male[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_geo_AUS_male_ncomp_6 = alpha_para_JS_div_geo_AUS_male_value_ncomp_6 = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select, interval = c(0, 0.5), dat = male_pop, horizon = iwk, criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_male_ncomp_6[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_value_ncomp_6[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

