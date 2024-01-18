###############################
# Numbers of factors selection
###############################

### The method of den Reijer et al. 2021 ###

k_criterion <- function(eigenvalue, plot = TRUE)
{
    n = length(eigenvalue)
    Hn = sum(1/(1:n))
    
    eigen_diff = eigenvalue[-n] - eigenvalue[-1]
    lambda_bar = 1/((2:n)*Hn)
    
    n_cross = which(eigen_diff <= lambda_bar)[1]
    
    if(n_cross > 1)
    {
        n_return = n_cross - 1
    } 
    else 
    {
        n_return = 1
    }
    
    if(plot)
    {
        plot(eigen_diff[1:(n_cross+6)], type = "b", xlab = expression("Eigenvalue number" ~ italic(k)), ylab = expression("Eigenvalue" ~~~ lambda))
        lines(lambda_bar[1:(n_cross+6)], type ="b", col = 2)
        abline(v = n_return, lty = 2, col = 4)
        legend("topright", lty = c(1,1), lwd = c(2,2), col = c(1,2), c(expression(lambda ~~ "(eigenvalue)  "), expression(bar(lambda) ~~ "(threshold)  ")))
    }
    return(n_return)
}

### Eigenratio $k$ selection method ###

select_K <- function(tau, eigenvalue)
{
    k_max = length(eigenvalue)
    k_all = rep(0, k_max-1)
    for(k in 1:(k_max-1))
    {
        k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
    }
    K_hat = which.min(k_all)
    return(K_hat)
}

