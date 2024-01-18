###########
# testing
###########

alpha_testing_h20 <- function(alpha_para, dat, method_ncomp, horizon, criterion, mat_Helmert)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 21 - horizon)
    for(iw in 1:(21 - horizon))
    {
        den_fore[,iw] = alpha_fun(object = dat[1:(n - 21 + iw),], alpha_val = alpha_para, ncomp_method = method_ncomp, 
                                  fh = horizon, Helmert_mat = mat_Helmert)[,horizon]
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

#########
# KL_div
#########

## female

# eigenvalue

alpha_KL_div_female_testing_h20 = clr_KL_div_female_testing_h20 =
ilr_KL_div_female_testing_h20 = eda_KL_div_female_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_KL_div_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_female_h20[iwk], 
                                                         dat = female_pop, method_ncomp = "eigenvalue", 
                                                         horizon = iwk, criterion = "KL_div", 
                                                         mat_Helmert = TRUE)
    
    ilr_KL_div_female_testing_h20[iwk]  = alpha_testing_h20(alpha_para = 0, dat = female_pop, horizon = iwk, 
                                                        criterion = "KL_div", method_ncomp = "eigenvalue", 
                                                        mat_Helmert = TRUE)
        
    clr_KL_div_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, horizon = iwk,
                                                           criterion = "KL_div", method_ncomp = "eigenvalue",
                                                           mat_Helmert = FALSE)
    
    eda_KL_div_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = 1, dat = female_pop, horizon = iwk,
                                                       criterion = "KL_div", method_ncomp = "eigenvalue", 
                                                       mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed (K = 6)

alpha_KL_div_female_testing_h20_ncomp_6 = clr_KL_div_female_testing_h20_ncomp_6 =
ilr_KL_div_female_testing_h20_ncomp_6 = eda_KL_div_female_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_KL_div_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_female_ncomp_6_h20[iwk], 
                                                                 dat = female_pop, method_ncomp = "fixed", 
                                                                 horizon = iwk, criterion = "KL_div",
                                                                 mat_Helmert = TRUE)
    
    ilr_KL_div_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, horizon = iwk,
                                                               criterion = "KL_div", method_ncomp = "fixed",
                                                               mat_Helmert = TRUE)

    clr_KL_div_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, horizon = iwk,
                                                                   criterion = "KL_div", method_ncomp = "fixed",
                                                                   mat_Helmert = FALSE)
    
    eda_KL_div_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = female_pop, horizon = iwk,
                                                               criterion = "KL_div", method_ncomp = "fixed",
                                                               mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

## male

# eigenvalue

alpha_KL_div_male_testing_h20 = clr_KL_div_male_testing_h20 =
ilr_KL_div_male_testing_h20 = eda_KL_div_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_KL_div_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_male_h20[iwk], 
                                                         dat = male_pop, method_ncomp = "eigenvalue", 
                                                         horizon = iwk, criterion = "KL_div", 
                                                         mat_Helmert = TRUE)
    
    ilr_KL_div_male_testing_h20[iwk]  = alpha_testing_h20(alpha_para = 0, dat = male_pop, horizon = iwk, 
                                                        criterion = "KL_div", method_ncomp = "eigenvalue", 
                                                        mat_Helmert = TRUE)
    
    clr_KL_div_male_testing_h20[iwk]  = alpha_testing_h20(alpha_para = 0, dat = male_pop, horizon = iwk,
                                                        criterion = "KL_div", method_ncomp = "eigenvalue",
                                                        mat_Helmert = FALSE)
    
    eda_KL_div_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = 1, dat = male_pop, horizon = iwk,
                                                       criterion = "KL_div", method_ncomp = "eigenvalue", 
                                                       mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed

alpha_KL_div_male_testing_h20_ncomp_6 = ilr_KL_div_male_testing_h20_ncomp_6 = eda_KL_div_male_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_KL_div_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_male_ncomp_6_h20[iwk], 
                                                                 dat = male_pop, method_ncomp = "fixed", 
                                                                 horizon = iwk, criterion = "KL_div",
                                                                 mat_Helmert = TRUE)
    
    ilr_KL_div_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = male_pop, horizon = iwk,
                                                               criterion = "KL_div", method_ncomp = "fixed",
                                                               mat_Helmert = TRUE)
    
    eda_KL_div_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = male_pop, horizon = iwk,
                                                               criterion = "KL_div", method_ncomp = "fixed",
                                                               mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

KL_div_testing_h20_results = cbind(alpha_KL_div_female_testing_h20, ilr_KL_div_female_testing_h20, eda_KL_div_female_testing_h20,
                               alpha_KL_div_male_testing_h20,   ilr_KL_div_male_testing_h20,   eda_KL_div_male_testing_h20)

KL_div_testing_h20_results_ncomp_6 = cbind(alpha_KL_div_female_testing_h20_ncomp_6, ilr_KL_div_female_testing_h20_ncomp_6, eda_KL_div_female_testing_h20_ncomp_6, 
                                       alpha_KL_div_male_testing_h20_ncomp_6,   ilr_KL_div_male_testing_h20_ncomp_6,   eda_KL_div_male_testing_h20_ncomp_6)    

colnames(KL_div_testing_h20_results) = colnames(KL_div_testing_h20_results_ncomp_6) = c("alpha (F)", "ilr (F)", "eda (F)",
                                                                                        "alpha (M)", "ilr (M)", "eda (M)")

round(colMeans(KL_div_testing_h20_results), 4)
round(colMeans(KL_div_testing_h20_results_ncomp_6), 4)

################
# JS_div_simple
################

## female

# eigenvalue

alpha_JS_div_simple_female_testing_h20 = ilr_JS_div_simple_female_testing_h20 = eda_JS_div_simple_female_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_simple_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_simple_AUS_female_h20[iwk], 
                                                            dat = female_pop, method_ncomp = "eigenvalue", 
                                                            horizon = iwk, criterion = "JS_div_simple", 
                                                            mat_Helmert = TRUE)
    
    ilr_JS_div_simple_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, 
                                                          method_ncomp = "eigenvalue", horizon = iwk, 
                                                          criterion = "JS_div_simple", mat_Helmert = TRUE)
    
    eda_JS_div_simple_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = 1, dat = female_pop, 
                                                          method_ncomp = "eigenvalue", horizon = iwk, 
                                                          criterion = "JS_div_simple", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed

alpha_JS_div_simple_female_testing_h20_ncomp_6 = ilr_JS_div_simple_female_testing_h20_ncomp_6 = eda_JS_div_simple_female_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_simple_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_simple_AUS_female_ncomp_6_h20[iwk],
                                                                    dat = female_pop, method_ncomp = "fixed", horizon = iwk,
                                                                    criterion = "JS_div_simple", mat_Helmert = TRUE)
    
    ilr_JS_div_simple_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, 
                                                                  method_ncomp = "fixed", horizon = iwk,
                                                                  criterion = "JS_div_simple", mat_Helmert = TRUE)
    
    eda_JS_div_simple_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = female_pop, 
                                                                  method_ncomp = "fixed", horizon = iwk,
                                                                  criterion = "JS_div_simple", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

## male

# eigenvalue

alpha_JS_div_simple_male_testing_h20 = clr_JS_div_simple_male_testing_h20 =
ilr_JS_div_simple_male_testing_h20 = eda_JS_div_simple_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    # alpha_para_JS_div_simple_AUS_male_h20
    
    alpha_JS_div_simple_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_male_h20[iwk],
                                                                dat = male_pop, method_ncomp = "eigenvalue",
                                                                horizon = iwk, criterion = "JS_div_simple", 
                                                                mat_Helmert = TRUE)
    
    ilr_JS_div_simple_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = 0, dat = male_pop, 
                                                              method_ncomp = "eigenvalue", horizon = iwk, 
                                                              criterion = "JS_div_simple", mat_Helmert = TRUE)
                                                              
    clr_JS_div_simple_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = 0, dat = male_pop,
                                                                method_ncomp = "eigenvalue", horizon = iwk,
                                                                criterion = "JS_div_simple", mat_Helmert = FALSE)
    
    eda_JS_div_simple_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = 1, dat = male_pop, 
                                                              method_ncomp = "eigenvalue", horizon = iwk, 
                                                              criterion = "JS_div_simple", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed

alpha_JS_div_simple_male_testing_h20_ncomp_6 = ilr_JS_div_simple_male_testing_h20_ncomp_6 = eda_JS_div_simple_male_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_simple_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_simple_AUS_male_ncomp_6_h20[iwk],
                                                                        dat = male_pop, method_ncomp = "fixed", horizon = iwk,
                                                                        criterion = "JS_div_simple", mat_Helmert = TRUE)
    
    ilr_JS_div_simple_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = male_pop, 
                                                                      method_ncomp = "fixed", horizon = iwk,
                                                                      criterion = "JS_div_simple", mat_Helmert = TRUE)
    
    eda_JS_div_simple_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = male_pop, 
                                                                      method_ncomp = "fixed", horizon = iwk,
                                                                      criterion = "JS_div_simple", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

JS_div_simple_testing_h20_results = cbind(alpha_JS_div_simple_female_testing_h20, ilr_JS_div_simple_female_testing_h20, eda_JS_div_simple_female_testing_h20,
                                          alpha_JS_div_simple_male_testing_h20,   ilr_JS_div_simple_male_testing_h20,   eda_JS_div_simple_male_testing_h20)

JS_div_simple_testing_h20_results_ncomp_6 = cbind(alpha_JS_div_simple_female_testing_h20_ncomp_6, ilr_JS_div_simple_female_testing_h20_ncomp_6, 		
                          											  eda_JS_div_simple_female_testing_h20_ncomp_6, alpha_JS_div_simple_male_testing_h20_ncomp_6,   
                          											  ilr_JS_div_simple_male_testing_h20_ncomp_6,   eda_JS_div_simple_male_testing_h20_ncomp_6)

colnames(JS_div_simple_testing_h20_results) = colnames(JS_div_simple_testing_h20_results_ncomp_6) = c("alpha (F)", "ilr (F)", "eda (F)", 
                                                                                                      "alpha (M)", "ilr (M)", "eda (M)")

round(colMeans(JS_div_simple_testing_h20_results), 4)
round(colMeans(JS_div_simple_testing_h20_results_ncomp_6), 4)

#############
# JS_div_geo
#############

## female

# eigenvalue

alpha_JS_div_geo_female_testing_h20 = ilr_JS_div_geo_female_testing_h20 = eda_JS_div_geo_female_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_geo_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_geo_AUS_female_h20[iwk], 
                                                             dat = female_pop, method_ncomp = "eigenvalue", 
                                                             horizon = iwk, criterion = "JS_div_geo",
                                                             mat_Helmert = TRUE)
    
    ilr_JS_div_geo_female_testing_h20[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, 
                                                           method_ncomp = "eigenvalue", horizon = iwk, criterion = "JS_div_geo",
                                                           mat_Helmert = TRUE)
        
    eda_JS_div_geo_female_testing_h20[iwk]   = alpha_testing_h20(alpha_para = 1, dat = female_pop, 
                                                             method_ncomp = "eigenvalue", horizon = iwk, 
                                                             criterion = "JS_div_geo", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed

alpha_JS_div_geo_female_testing_h20_ncomp_6 = ilr_JS_div_geo_female_testing_h20_ncomp_6 = eda_JS_div_geo_female_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_geo_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_geo_AUS_female_ncomp_6_h20[iwk],
                                                                     dat = female_pop, method_ncomp = "fixed", horizon = iwk,
                                                                     criterion = "JS_div_geo", mat_Helmert = TRUE)
    
    ilr_JS_div_geo_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = female_pop, 
                                                                   method_ncomp = "fixed", horizon = iwk,
                                                                   criterion = "JS_div_geo", mat_Helmert = TRUE)
        
    eda_JS_div_geo_female_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = female_pop,
                                                                   method_ncomp = "fixed", horizon = iwk,
                                                                   criterion = "JS_div_geo", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

## male

# eigenvalue

alpha_JS_div_geo_male_testing_h20 = ilr_JS_div_geo_male_testing_h20 = eda_JS_div_geo_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    # alpha_para_JS_div_geo_AUS_male_h20
    alpha_JS_div_geo_male_testing_h20[iwk] = alpha_testing_h20(alpha_para = alpha_para_KL_div_AUS_male_h20[iwk],
                                                           dat = male_pop, method_ncomp = "eigenvalue",
                                                           horizon = iwk, criterion = "JS_div_geo", mat_Helmert = TRUE)
    
    ilr_JS_div_geo_male_testing_h20[iwk]  = alpha_testing_h20(alpha_para = 0, dat = male_pop,   
                                                          method_ncomp = "eigenvalue", horizon = iwk, 
                                                          criterion = "JS_div_geo", mat_Helmert = TRUE)
    
    eda_JS_div_geo_male_testing_h20[iwk]   = alpha_testing_h20(alpha_para = 1, dat = male_pop,  
                                                           method_ncomp = "eigenvalue", horizon = iwk, 
                                                           criterion = "JS_div_geo", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}

# fixed

alpha_JS_div_geo_male_testing_h20_ncomp_6 = ilr_JS_div_geo_male_testing_h20_ncomp_6 = eda_JS_div_geo_male_testing_h20_ncomp_6 = vector("numeric", 20)
for(iwk in 1:20)
{
    alpha_JS_div_geo_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = alpha_para_JS_div_geo_AUS_male_ncomp_6_h20[iwk],
                                                               dat = male_pop, method_ncomp = "fixed", horizon = iwk,
                                                               criterion = "JS_div_geo", mat_Helmert = TRUE)
    
    ilr_JS_div_geo_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 0, dat = male_pop,
                                                             method_ncomp = "fixed", horizon = iwk,
                                                             criterion = "JS_div_geo", mat_Helmert = TRUE)
    
    eda_JS_div_geo_male_testing_h20_ncomp_6[iwk] = alpha_testing_h20(alpha_para = 1, dat = male_pop,
                                                             method_ncomp = "fixed", horizon = iwk,
                                                             criterion = "JS_div_geo", mat_Helmert = TRUE)
    print(iwk); rm(iwk)
}


JS_div_geo_testing_h20_results = cbind(alpha_JS_div_geo_female_testing_h20, ilr_JS_div_geo_female_testing_h20, eda_JS_div_geo_female_testing_h20,
                                       alpha_JS_div_geo_male_testing_h20,   ilr_JS_div_geo_male_testing_h20,   eda_JS_div_geo_male_testing_h20)

JS_div_geo_testing_h20_results_ncomp_6 = cbind(alpha_JS_div_geo_female_testing_h20_ncomp_6, ilr_JS_div_geo_female_testing_h20_ncomp_6, eda_JS_div_geo_female_testing_h20_ncomp_6, 
                                               alpha_JS_div_geo_male_testing_h20_ncomp_6,   ilr_JS_div_geo_male_testing_h20_ncomp_6,   eda_JS_div_geo_male_testing_h20_ncomp_6)

colnames(JS_div_geo_testing_h20_results) = colnames(JS_div_geo_testing_h20_results_ncomp_6) = c("alpha (F)", "clr (F)", "eda (F)", 
                                                                                                "alpha (M)", "clr (M)", "eda (M)")

xtable(rbind(JS_div_geo_testing_h20_results, colMeans(JS_div_geo_testing_h20_results)), digits = 4)
xtable(rbind(JS_div_geo_testing_h20_results_ncomp_6, colMeans(JS_div_geo_testing_h20_results_ncomp_6)), digits = 4)
