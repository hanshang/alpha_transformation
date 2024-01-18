################
# Package names
################

#devtools::install_github("mpascariu/MortalityForecast")

packages <- c("Compositional", "psych", "ftsa", "meboot", "pracma", "reldist", "flexmix", "demography", 
              "MortalityForecast", "MortalityLaws", "xtable")

# Install packages not yet installed

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading

invisible(lapply(packages, library, character.only = TRUE))

##########################
# set a working directory
##########################

age = 0:110
year = 1921:2020
n_year = length(year)
female_qx = t(matrix(read.table("AUS_lt_female_death.txt", header = TRUE)[,"qx"], 111, n_year))
male_qx   = t(matrix(read.table("AUS_lt_male_death.txt", header = TRUE)[,"qx"], 111, n_year))

n_col = ncol(female_qx)
n_row = nrow(female_qx)
female_pop = male_pop = matrix(NA, n_row, n_col)
for(ij in 1:n_row)
{
    start_pop_female = start_pop_male = 10^5
    for(ik in 1:n_col)
    {
        female_pop[ij,ik] = female_qx[ij,ik] * start_pop_female
        start_pop_female  = start_pop_female - female_pop[ij,ik]
        
        male_pop[ij,ik] = male_qx[ij,ik] * start_pop_male
        start_pop_male  = start_pop_male - male_pop[ij,ik]
        rm(ik)
    }
    print(ij); rm(ij)
}
rownames(female_pop) = rownames(male_pop) = year
colnames(female_pop) = colnames(male_pop) = age

# plot figures

par(mfrow = c(1, 2))
plot(fts(age, t(female_pop)), xlab = "Age", ylab = "Life-table death counts", 
     main = "Australian female data\n (1921-2020)", ylim = c(0, 0.08))
plot(fts(age, t(male_pop)),   xlab = "Age", ylab = "", main = "Australian male data\n (1921-2020)",
     ylim = c(0, 0.08))

# save figures

plot(fts(age, t(female_pop)), xlab = "Age", ylab = "Life-table death counts",
     main = "Australia: female data\n (1921-2020)", ylim = c(0, 0.08 * 10^5))
plot(fts(age, t(male_pop)),   xlab = "Age", ylab = "", main = "Australia: male data\n (1921-2020)",
     ylim = c(0, 0.08* 10^5))
