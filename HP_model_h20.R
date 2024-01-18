###########################
# Point forecasts (h = 20)
###########################

HP_testing_h20 <- function(dat, horizon, criterion, no_moment, PI_level)
{
    n = nrow(dat)
    age = 0:(ncol(dat)-1)

    mortlaw_dat_fore = matrix(NA, ncol(dat), 21 - horizon)
    for(iw in 1:(21 - horizon))
    {
        model_MEM = model.MEM(data = t(dat[1:(n - 21 + iw),]), x = age, y = (1:n)[1:(n - 21 + iw)], n = no_moment)
        mortlaw_dat_fore[,iw] = predict(model_MEM, h = 20, level = PI_level)$predicted.values[,horizon]
        rm(iw); rm(model_MEM)
    }

    err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        data_c = cbind(True = dat[(n - 21 + horizon + ik),], forecast = as.numeric(mortlaw_dat_fore[,ik]))
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
                colnames(P_M) = colnames(E_M) = c("True", "M")
                err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
                rm(M); rm(P_M); rm(E_M)
            }
            else if(criterion == "JS_div_geo")
            {
                M = apply(data_c, 1, geometric.mean)
                P_M = cbind(data_c[,1], M)
                E_M = cbind(as.numeric(data_c[,2]), M)
                colnames(P_M) = colnames(E_M) = c("True", "M")
                err[ik] = as.numeric(0.5 * KLdiv(P_M) + 0.5 * KLdiv(E_M))[3]
                rm(M); rm(P_M); rm(E_M)
            }
            else
            {
                warning("Please specify a criterion from the list.")
            }
        }
    }
    rm(n); rm(mortlaw_dat_fore)
    return(mean(as.numeric(err), na.rm = TRUE))
}

# point forecasts

HP_KL_div_female_testing_h20 = HP_KL_div_male_testing_h20 =
HP_JS_div_simple_female_testing_h20 = HP_JS_div_simple_male_testing_h20 =
HP_JS_div_geo_female_testing_h20 = HP_JS_div_geo_male_testing_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    # KL_div

    HP_KL_div_female_testing_h20[iwk] = HP_testing_h20(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                                       no_moment = 4, PI_level = 95)

    HP_KL_div_male_testing_h20[iwk]   = HP_testing_h20(dat = male_pop,   horizon = iwk, criterion = "KL_div",
                                                       no_moment = 4, PI_level = 95)

    # JS_div_simple

    HP_JS_div_simple_female_testing_h20[iwk] = HP_testing_h20(dat = female_pop, horizon = iwk,
                                                              criterion = "JS_div_simple", no_moment = 4,
                                                              PI_level = 95)

    HP_JS_div_simple_male_testing_h20[iwk]   = HP_testing_h20(dat = male_pop,   horizon = iwk,
                                                              criterion = "JS_div_simple", no_moment = 4,
                                                              PI_level = 95)

    # JS_div_geo

    HP_JS_div_geo_female_testing_h20[iwk] = HP_testing_h20(dat = female_pop, horizon = iwk,
                                                           criterion = "JS_div_geo", no_moment = 4,
                                                           PI_level = 95)

    HP_JS_div_geo_male_testing_h20[iwk]   = HP_testing_h20(dat = male_pop,   horizon = iwk,
                                                           criterion = "JS_div_geo", no_moment = 4,
                                                           PI_level = 95)
    print(iwk); rm(iwk)
}

HP_results_testing_h20 = cbind(HP_KL_div_female_testing_h20, HP_JS_div_simple_female_testing_h20, HP_JS_div_geo_female_testing_h20,
                               HP_KL_div_male_testing_h20,   HP_JS_div_simple_male_testing_h20,   HP_JS_div_geo_male_testing_h20)
colnames(HP_results_testing_h20) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                     "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(HP_results_testing_h20) = 1:20

#############################
# interval forecasts (h = 20)
#############################

HP_testing_interval_h20 <- function(dat, horizon, criterion, no_moment, PI_level)
{
    n = nrow(dat)
    age = 0:(ncol(dat) - 1)
    
    mortlaw_dat_fore_lb = mortlaw_dat_fore_ub = matrix(NA, ncol(dat), 21 - horizon)
    for(iw in 1:(21 - horizon))
    {
        model_MEM = model.MEM(data = t(dat[1:(n - 21 + iw),]), x = age, y = (1:n)[1:(n - 21 + iw)], n = no_moment)
        predict_model_MEM = predict(model_MEM, h = 20, level = PI_level)$conf.intervals$predicted.values
        if(PI_level == 95)
        {
            mortlaw_dat_fore_lb[,iw] = predict_model_MEM$L95[,horizon]
            mortlaw_dat_fore_ub[,iw] = predict_model_MEM$U95[,horizon]
        }
        else if(PI_level == 80)
        {
            mortlaw_dat_fore_lb[,iw] = predict_model_MEM$L80[,horizon]
            mortlaw_dat_fore_ub[,iw] = predict_model_MEM$U80[,horizon]
        }
        else
        {
            warning("PI level should either be 95 or 80.")
        }
        rm(iw); rm(model_MEM)
    }
    
    int_err = vector("numeric", (21 - horizon))
    for(ik in 1:(21 - horizon))
    {
        if(criterion == "CPD")
        {
            int_err[ik] = cpd(holdout = dat[(n - 21 + horizon + ik),], lb = mortlaw_dat_fore_lb[,ik],
                              ub = mortlaw_dat_fore_ub[,ik], alpha = (100 - PI_level)/100)
        }
        else if(criterion == "score")
        {
            int_err[ik] = interval_score(holdout = dat[(n - 21 + horizon + ik),], lb = mortlaw_dat_fore_lb[,ik],
                                         ub = mortlaw_dat_fore_ub[,ik], alpha = (100 - PI_level)/100)
        }
        else
        {
            warning("Criterion must be either CPD or interval score")
        }
        rm(ik)
    }
    result_ave = mean(int_err)
    rm(int_err); rm(n)
    return(result_ave)
}

## female

# CPD

HP_testing_cpd_PI_80_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_cpd_PI_80_h20[iwk] = HP_testing_interval_h20(dat = female_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_cpd_PI_95_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_cpd_PI_95_h20[iwk] = HP_testing_interval_h20(dat = female_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

# interval score

HP_testing_score_PI_80_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_score_PI_80_h20[iwk] = HP_testing_interval_h20(dat = female_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_score_PI_95_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_score_PI_95_h20[iwk] = HP_testing_interval_h20(dat = female_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

HP_testing_female_h20 = cbind(HP_testing_cpd_PI_80_h20, HP_testing_score_PI_80_h20, HP_testing_cpd_PI_95_h20, HP_testing_score_PI_95_h20)
colnames(HP_testing_female_h20) = c("CPD (80)", "score (80)", "CPD (95)", "score (95)")

round(colMeans(HP_testing_female_h20), 4) # 0.8199 0.0149 0.8513 0.0828

## male

# CPD

HP_testing_male_cpd_PI_80_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_male_cpd_PI_80_h20[iwk] = HP_testing_interval_h20(dat = male_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_male_cpd_PI_95_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_male_cpd_PI_95_h20[iwk] = HP_testing_interval_h20(dat = male_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 95)
}

# interval score

HP_testing_male_score_PI_80_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_male_score_PI_80_h20[iwk] = HP_testing_interval_h20(dat = male_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_male_score_PI_95_h20 = vector("numeric", 20)
for(iwk in 1:20)
{
    HP_testing_male_score_PI_95_h20[iwk] = HP_testing_interval_h20(dat = male_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

HP_testing_male_h20 = cbind(HP_testing_male_cpd_PI_80_h20, HP_testing_male_score_PI_80_h20, HP_testing_male_cpd_PI_95_h20, HP_testing_male_score_PI_95_h20)
colnames(HP_testing_male_h20) = c("CPD (80)", "score (80)", "CPD (95)", "score (95)")

round(colMeans(HP_testing_male_h20), 4) # 0.8443 0.0144 0.9239 0.0794

