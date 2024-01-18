#########################
# MOM (4 finite moments)
#########################

HP_testing <- function(dat, horizon, criterion, no_moment, PI_level)
{
    n = nrow(dat)
    age = 0:(ncol(dat)-1)

    mortlaw_dat_fore = matrix(NA, ncol(dat), 11 - horizon)
    for(iw in 1:(11 - horizon))
    {
        model_MEM = model.MEM(data = t(dat[1:(n - 11 + iw),]), x = age, y = (1:n)[1:(n - 11 + iw)], n = no_moment)
        mortlaw_dat_fore[,iw] = predict(model_MEM, h = 10, level = PI_level)$predicted.values[,horizon]
        rm(iw); rm(model_MEM)
    }

    err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        data_c = cbind(True = dat[(n - 11 + horizon + ik),], forecast = as.numeric(mortlaw_dat_fore[,ik]))
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

HP_KL_div_female_testing = HP_KL_div_male_testing =
HP_JS_div_simple_female_testing = HP_JS_div_simple_male_testing =
HP_JS_div_geo_female_testing = HP_JS_div_geo_male_testing = vector("numeric", 10)
for(iwk in 1:10)
{
    # KL_div

    HP_KL_div_female_testing[iwk] = HP_testing(dat = female_pop, horizon = iwk, criterion = "KL_div",
                                               no_moment = 4, PI_level = 95)

    HP_KL_div_male_testing[iwk]   = HP_testing(dat = male_pop,   horizon = iwk, criterion = "KL_div",
                                               no_moment = 4, PI_level = 95)

    # JS_div_simple

    HP_JS_div_simple_female_testing[iwk] = HP_testing(dat = female_pop, horizon = iwk,
                                                      criterion = "JS_div_simple", no_moment = 4,
                                                      PI_level = 95)

    HP_JS_div_simple_male_testing[iwk]   = HP_testing(dat = male_pop,   horizon = iwk,
                                                      criterion = "JS_div_simple", no_moment = 4,
                                                      PI_level = 95)

    # JS_div_geo

    HP_JS_div_geo_female_testing[iwk] = HP_testing(dat = female_pop, horizon = iwk,
                                                   criterion = "JS_div_geo", no_moment = 4,
                                                   PI_level = 95)

    HP_JS_div_geo_male_testing[iwk]   = HP_testing(dat = male_pop,   horizon = iwk,
                                                   criterion = "JS_div_geo", no_moment = 4,
                                                   PI_level = 95)
    print(iwk); rm(iwk)
}

HP_results_testing = cbind(HP_KL_div_female_testing, HP_JS_div_simple_female_testing, HP_JS_div_geo_female_testing,
                           HP_KL_div_male_testing,   HP_JS_div_simple_male_testing,   HP_JS_div_geo_male_testing)
colnames(HP_results_testing) = c("KL_div (F)", "JS_div_simple (F)", "JS_div_geo (F)",
                                 "KL_div (M)", "JS_div_simple (M)", "JS_div_geo (M)")
rownames(HP_results_testing) = 1:10

#############################
# interval forecasts (h = 10)
#############################

HP_testing_interval <- function(dat, horizon, criterion, no_moment, PI_level)
{
    n = nrow(dat)
    age = 0:(ncol(dat) - 1)

    mortlaw_dat_fore_lb = mortlaw_dat_fore_ub = matrix(NA, ncol(dat), 11 - horizon)
    for(iw in 1:(11 - horizon))
    {
        model_MEM = model.MEM(data = t(dat[1:(n - 11 + iw),]), x = age, y = (1:n)[1:(n - 11 + iw)], n = no_moment)
        predict_model_MEM = predict(model_MEM, h = 10, level = PI_level)$conf.intervals$predicted.values
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

    int_err = vector("numeric", (11 - horizon))
    for(ik in 1:(11 - horizon))
    {
        if(criterion == "CPD")
        {
            int_err[ik] = cpd(holdout = dat[(n - 11 + horizon + ik),], lb = mortlaw_dat_fore_lb[,ik],
                              ub = mortlaw_dat_fore_ub[,ik], alpha = (100 - PI_level)/100)
        }
        else if(criterion == "score")
        {
            int_err[ik] = interval_score(holdout = dat[(n - 11 + horizon + ik),], lb = mortlaw_dat_fore_lb[,ik],
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

HP_testing_cpd_PI_80 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_cpd_PI_80[iwk] = HP_testing_interval(dat = female_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_cpd_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_cpd_PI_95[iwk] = HP_testing_interval(dat = female_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

# interval score

HP_testing_score_PI_80 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_score_PI_80[iwk] = HP_testing_interval(dat = female_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_score_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_score_PI_95[iwk] = HP_testing_interval(dat = female_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

HP_testing_female = cbind(HP_testing_cpd_PI_80, HP_testing_score_PI_80, HP_testing_cpd_PI_95, HP_testing_score_PI_95)
colnames(HP_testing_female) = c("CPD (80)", "score (80)", "CPD (95)", "score (95)")

round(colMeans(HP_testing_female), 4) # 0.8171 0.0100 0.8702 0.0543

## male

# CPD

HP_testing_male_cpd_PI_80 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_male_cpd_PI_80[iwk] = HP_testing_interval(dat = male_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_male_cpd_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_male_cpd_PI_95[iwk] = HP_testing_interval(dat = male_pop, horizon = iwk, criterion = "CPD", no_moment = 4, PI_level = 95)
}

# interval score

HP_testing_male_score_PI_80 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_male_score_PI_80[iwk] = HP_testing_interval(dat = male_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 80)
    print(iwk); rm(iwk)
}

HP_testing_male_score_PI_95 = vector("numeric", 10)
for(iwk in 1:10)
{
    HP_testing_male_score_PI_95[iwk] = HP_testing_interval(dat = male_pop, horizon = iwk, criterion = "score", no_moment = 4, PI_level = 95)
    print(iwk); rm(iwk)
}

HP_testing_male = cbind(HP_testing_male_cpd_PI_80, HP_testing_male_score_PI_80, HP_testing_male_cpd_PI_95, HP_testing_male_score_PI_95)
colnames(HP_testing_male) = c("CPD (80)", "score (80)", "CPD (95)", "score (95)")

round(colMeans(HP_testing_male), 4) # 0.8351 0.0096 0.9145 0.0512

