########
# plots
########

alpha_fun_ex <- function(object, alpha_val, ncomp_tuning = 0.001, fh, Helmert_mat)
{
    object_alpha = alfa(object, a = alpha_val, h = Helmert_mat)$aff
    SVD_decomp = svd(object_alpha)
    ncomp = select_K(tau = ncomp_tuning, SVD_decomp$d^2)
    
    basis = as.matrix(SVD_decomp$v[,1:ncomp])
    score = object_alpha %*% basis
    recon = t(basis %*% t(score))
    
    # reconstruction
    
    object_recon = alfainv(recon, a = alpha_val, h = Helmert_mat)
    R2 = round(1 - sum((object_recon - object)^2)/sum((object - colMeans(object))^2), 4)
    RMSE = ftsa:::rmse(forecast = object_recon, true = object)
      
    # forecasts of principal component scores
    
    score_fore = matrix(NA, ncomp, fh)
    for(ik in 1:ncomp)
    {
        score_fore[ik,] = forecast(auto.arima(as.numeric(score[,ik])), h = fh)$mean
    }
    
    # obtain forecasts in real-valued space
    
    fore_val = t(basis %*% score_fore)
    object_fore = alfainv(fore_val, a = alpha_val)
    return(list(forecast_val = t(object_fore), object_alpha = object_alpha, ncomp = ncomp, 
                basis = basis, score = score, score_fore = score_fore, R2 = R2, RMSE = RMSE))
}

female_pop_ex = alpha_fun_ex(object = female_pop, alpha_val = alpha_para_KL_div_AUS_female[10], 
                             fh = 10, Helmert_mat = TRUE)
male_pop_ex   = alpha_fun_ex(object = male_pop,   alpha_val = alpha_para_KL_div_AUS_male[10], 
                             fh = 10, Helmert_mat = TRUE)
colnames(female_pop_ex$forecast_val) = colnames(male_pop_ex$forecast_val) = 2021:2030

female_pop_ex$R2  # 0.9968
male_pop_ex$R2    # 0.9915

round(female_pop_ex$RMSE, 5) # 0.00091
round(male_pop_ex$RMSE, 5)   # 0.00138    

female_pop_ex$ncomp # 2
male_pop_ex$ncomp   # 2 

# Helmert 0 and 1

female_pop_ex_bench_zero = alpha_fun_ex(object = female_pop, alpha_val = 0, fh = 10, Helmert_mat = TRUE)
female_pop_ex_bench_one  = alpha_fun_ex(object = female_pop, alpha_val = 1, fh = 10, Helmert_mat = TRUE)
  
plot(fts(head(age, -1), t(female_pop_ex_bench_zero$object_alpha), xname = "Age", yname = expression(bold(z))),
      main = "Australia: female data (1921-2020)")

plot(fts(head(age, -1), t(female_pop_ex_bench_one$object_alpha), xname = "Age", yname = expression(bold(z))),
     main = "Australia: female data (1921-2020)")

male_pop_ex_bench_zero = alpha_fun_ex(object = male_pop, alpha_val = 0, fh = 10, Helmert_mat = TRUE)
male_pop_ex_bench_one  = alpha_fun_ex(object = male_pop, alpha_val = 1, fh = 10, Helmert_mat = TRUE)

colnames(female_pop_ex_bench_zero$forecast_val) = colnames(male_pop_ex_bench_zero$forecast_val) = 
colnames(female_pop_ex_bench_one$forecast_val) = colnames(male_pop_ex_bench_one$forecast_val) = 2021:2030

female_pop_ex_bench_zero$R2 # 0.9953
male_pop_ex_bench_zero$R2   # 0.9911

round(female_pop_ex_bench_zero$RMSE, 5) # 0.00111
round(male_pop_ex_bench_zero$RMSE, 5)   # 0.00141

female_pop_ex_bench_zero$ncomp # 2
male_pop_ex_bench_zero$ncomp   # 2 

####################
## plotting figures
####################

# female

plot(fts(head(age, -1), t(female_pop_ex$object_alpha), xname = "Age", yname = expression(bold(z))),
     main = "Australia: female data (1921-2020)")
plot(head(age, -1), female_pop_ex$basis[,1], type = "l", xlab = "Age", ylab = expression(hat(phi)[1]))
plot(head(age, -1), female_pop_ex$basis[,2], type = "l", xlab = "Age", ylab = expression(hat(phi)[2]))

plot(forecast(auto.arima(ts(as.numeric(female_pop_ex$score[,1]), start = head(year, 1), end = tail(year, 1))), h = 10),
     xlab = "Year", ylab = expression(hat(beta)[1]))
plot(forecast(auto.arima(ts(as.numeric(female_pop_ex$score[,2]), start = head(year, 1), end = tail(year, 1))), h = 10),
     xlab = "Year", ylab = expression(hat(beta)[2]))
plot(fts(age, female_pop_ex$forecast_val, xname = "Age", yname = "Life-table death count"),
     main = "Australia: female data (2021-2030)")

# male

plot(fts(head(age, -1), t(male_pop_ex$object_alpha), xname = "Age", yname = expression(bold(z))),
     main = "Australia: male data (1921-2020)")
plot(head(age, -1), male_pop_ex$basis[,1], type = "l", xlab = "Age", ylab = expression(hat(phi)[1]))
plot(head(age, -1), male_pop_ex$basis[,2], type = "l", xlab = "Age", ylab = expression(hat(phi)[2]))

plot(forecast(auto.arima(ts(as.numeric(male_pop_ex$score[,1]), start = head(year, 1), end = tail(year, 1))), h = 10),
     xlab = "Year", ylab = expression(hat(beta)[1]))
plot(forecast(auto.arima(ts(as.numeric(male_pop_ex$score[,2]), start = head(year, 1), end = tail(year, 1))), h = 10),
     xlab = "Year", ylab = expression(hat(beta)[2]))
plot(fts(age, male_pop_ex$forecast_val, xname = "Age", yname = "Life-table death count"),
     main = "Australia: male data (2021-2030)")
