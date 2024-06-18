require(Compositional)
require(psych)
require(flexmix)

#################
# rolling window
#################

# alpha_para: alpha parameter value
# dat: life-table death counts
# horizon: forecast horizon
# criterion: point forecast evaluation accuracy
# mat_Helmert: TRUE or FALSE
# method_ncomp: method for selecting the number of retained components

alpha_para_select_fixed_window <- function(alpha_para, dat, horizon, criterion, mat_Helmert, method_ncomp)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 11 - horizon)
    for(iw in 1:(11 - horizon))
    {
        den_fore[,iw] = alpha_fun(object = dat[iw:(n - 21 + iw),], alpha_val = alpha_para,
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

alpha_para_KL_div_AUS_female_fixed_window = alpha_para_KL_div_AUS_female_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk,
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_female_fixed_window[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# ncomp = 6

alpha_para_KL_div_AUS_female_ncomp_6_fixed_window = alpha_para_KL_div_AUS_female_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk, 
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_female_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# AUS_male

alpha_para_KL_div_AUS_male_fixed_window = alpha_para_KL_div_AUS_male_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk,  
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_male_fixed_window[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# ncomp = 6

alpha_para_KL_div_AUS_male_ncomp_6_fixed_window = alpha_para_KL_div_AUS_male_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk, criterion = "KL_div", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_male_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

#################
## JS_div_simple
#################

## Helmert = TRUE (sum to 0)

# AUS_female

alpha_para_JS_div_simple_AUS_female_fixed_window = alpha_para_JS_div_simple_AUS_female_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_simple_AUS_female_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_simple_AUS_female_ncomp_6_fixed_window = alpha_para_JS_div_simple_AUS_female_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_simple_AUS_female_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# AUS_male

alpha_para_JS_div_simple_AUS_male_fixed_window = alpha_para_JS_div_simple_AUS_male_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk,  criterion = "JS_div_simple",
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_simple_AUS_male_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_simple_AUS_male_ncomp_6_fixed_window = alpha_para_JS_div_simple_AUS_male_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk, criterion = "JS_div_simple", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_simple_AUS_male_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

#############
# JS_div_geo
#############

## Helmert = TRUE 

# AUS_female

alpha_para_JS_div_geo_AUS_female_fixed_window = alpha_para_JS_div_geo_AUS_female_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk,  criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_female_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_geo_AUS_female_ncomp_6_fixed_window = alpha_para_JS_div_geo_AUS_female_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = female_pop, horizon = iwk, criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_female_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

## AUS_male

alpha_para_JS_div_geo_AUS_male_fixed_window = alpha_para_JS_div_geo_AUS_male_value_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk,  criterion = "JS_div_geo",
                   mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_male_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_value_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

alpha_para_JS_div_geo_AUS_male_ncomp_6_fixed_window = alpha_para_JS_div_geo_AUS_male_value_ncomp_6_fixed_window = vector("numeric", 10)
for(iwk in 1:10)
{
    dum = optimise(f = alpha_para_select_fixed_window, interval = c(0, 1), dat = male_pop, horizon = iwk, criterion = "JS_div_geo", 
                   mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_male_ncomp_6_fixed_window[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_value_ncomp_6_fixed_window[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

