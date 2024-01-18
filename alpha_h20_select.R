########################################################
# Training set 1:(n-41)
# Use validation set (n-40):(n-20) to evaluate accuracy
# Testing set (n-19):n to evaluate accuracy
########################################################

alpha_para_select_h20 <- function(alpha_para, dat, horizon, criterion, mat_Helmert, method_ncomp)
{
    n = nrow(dat)
    den_fore = matrix(NA, ncol(dat), 21 - horizon)
    for(iw in 1:(21 - horizon))
    {
      den_fore[,iw] = alpha_fun(object = dat[1:(n - 41 + iw),], alpha_val = alpha_para,
                                ncomp_method = method_ncomp,
                                fh = horizon, Helmert_mat = mat_Helmert)[,horizon]
    }
    err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
      data_c = cbind(True = dat[(n - 41 + horizon + ik),], forecast = as.numeric(den_fore[,ik]))
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

# eigenvalue ratio

alpha_para_KL_div_AUS_female_h20 = alpha_para_KL_div_AUS_female_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_female_h20[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed (K = 6)

alpha_para_KL_div_AUS_female_ncomp_6_h20 = alpha_para_KL_div_AUS_female_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_female_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_KL_div_AUS_female_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

## male

# eigenvalue ratio

alpha_para_KL_div_AUS_male_h20 = alpha_para_KL_div_AUS_male_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_KL_div_AUS_male_h20[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed

alpha_para_KL_div_AUS_male_ncomp_6_h20 = alpha_para_KL_div_AUS_male_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "KL_div", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_KL_div_AUS_male_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_KL_div_AUS_male_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

################
# JS_div_simple
################

## female

# eigenvalue

alpha_para_JS_div_simple_AUS_female_h20 = alpha_para_JS_div_simple_AUS_female_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "JS_div_simple", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_simple_AUS_female_h20[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed

alpha_para_JS_div_simple_AUS_female_ncomp_6_h20 = alpha_para_JS_div_simple_AUS_female_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "JS_div_simple", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_simple_AUS_female_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_female_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

## male

# eigenvalue

alpha_para_JS_div_simple_AUS_male_h20 = alpha_para_JS_div_simple_AUS_male_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{                 
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "JS_div_simple", mat_Helmert = TRUE, method_ncomp = "eigenvalue")             
    alpha_para_JS_div_simple_AUS_male_h20[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed

alpha_para_JS_div_simple_AUS_male_ncomp_6_h20 = alpha_para_JS_div_simple_AUS_male_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "JS_div_simple", method_ncomp = "fixed", mat_Helmert = TRUE)
    alpha_para_JS_div_simple_AUS_male_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_JS_div_simple_AUS_male_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

#############
# JS_div_geo
#############

## female

# eigenvalue

alpha_para_JS_div_geo_AUS_female_h20 = alpha_para_JS_div_geo_AUS_female_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "JS_div_geo", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_female_h20[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed

alpha_para_JS_div_geo_AUS_female_ncomp_6_h20 = alpha_para_JS_div_geo_AUS_female_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = female_pop, horizon = iwk, 
                   criterion = "JS_div_geo", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_female_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_female_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

## male

# eigenvalue

alpha_para_JS_div_geo_AUS_male_h20 = alpha_para_JS_div_geo_AUS_male_h20_value = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "JS_div_geo", mat_Helmert = TRUE, method_ncomp = "eigenvalue")
    alpha_para_JS_div_geo_AUS_male_h20[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_h20_value[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

# fixed

alpha_para_JS_div_geo_AUS_male_ncomp_6_h20 = alpha_para_JS_div_geo_AUS_male_value_ncomp_6_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    dum = optimise(f = alpha_para_select_h20, interval = c(0, 0.5), dat = male_pop, horizon = iwk, 
                   criterion = "JS_div_geo", mat_Helmert = TRUE, method_ncomp = "fixed")
    alpha_para_JS_div_geo_AUS_male_ncomp_6_h20[iwk] = dum$minimum
    alpha_para_JS_div_geo_AUS_male_value_ncomp_6_h20[iwk] = dum$objective
    print(iwk); rm(iwk); rm(dum)
}

