#--------------------------------------------------------------------------
# SET UP                          
#--------------------------------------------------------------------------
# Description -------------------------------------------------------------

# Analysis of 6-8 week baby checks

# Set working directory ---------------------------------------------------

setwd('K:/')

# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse, lubridate, 
               forestplot,
               epitools, RColorBrewer, pals,
               ggpubr, # ggarrange to combine graphs
               gghighlight, # highlighting specific groups in graphs
               lmtest, sandwich) # modified Poisson, sandwich and clustering

# Disable scientific notation ---------------------------------------------

options(scipen = 999)

# Write functions ----------------------------------------------------------

# create results table from modified poisson regression
extract_mpoisson_results <- function (x, x_sandwich, output_name, dataset, variable) {
  ethnicity <- names(coef(x))
  ethnicity <- sub("M_eth18_2011", "", ethnicity) # extracts all text after ')' from ethnicity to avoid having to recode later (excluding the intercept, the below lines of code are deal with this)
  ethnicity[1] <- levels(eval(substitute(variable), dataset))[1] # takes the reference level of a specified factor variable and replaces it with 'Intercept' in ethnicity
  estimate<- x_sandwich$RR
  lower <- x_sandwich$LCI
  upper <- x_sandwich$UCI
  p <- x_sandwich$P
  output <- cbind(ethnicity,estimate,lower, upper,p)
  output <- tibble::as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% dplyr::mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$rr_ci <- paste(output$estimate, output$ci, sep =",")
  output_table <- dplyr::select(output,ethnicity, rr_ci )
  assign(x = output_name, value = output, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
}

# calculate stratum specific CIs using sandwich estimators
get_sandwich_stratum_ci <- function (x_sandwich, ethnicity_est, ethnicity_region_int) {
  round(exp(confint(x_sandwich, level = 0.95)[ethnicity_est,] + confint(x_sandwich, level = 0.95)[ethnicity_region_int,]), 2)
}

#--------------------------------------------------------------------------
# DESCRIPTIVE ANALYSIS: PROVISION OF 6-8 WEEK BABY CHECK                        
#--------------------------------------------------------------------------

load(file="~/babycheck_mainanalysis.Rdata")
load(file="~/babycheck_sensanalysis.Rdata")

# Provision of baby check: main -------------------------------------------------

# main analysis, all checks
eth_denom_main_all <- nrow(babycheck_mainanalysis)
provision_main_all <- babycheck_mainanalysis %>% 
  group_by(babycheck) %>% 
  count() %>% 
  mutate(percent = round(n/eth_denom_main_all*100, 1)) %>% 
  mutate(prop = n/eth_denom_main_all) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/eth_denom_main_all))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/eth_denom_main_all))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(babycheck == 1) %>% 
  mutate(M_eth18_2011 = "All")
provision_main_all_table <- provision_main_all %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(`n(%, 95%CI)`, M_eth18_2011)  

eth_denom_main <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main <- babycheck_mainanalysis %>% 
  group_by(babycheck, M_eth18_2011) %>% 
  count() %>% 
  ungroup()
provision_main <- provision_main %>% 
  left_join(eth_denom_main, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(babycheck == 1)  
provision_main_table <- provision_main %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, babycheck)) 

provision_main_all <- rbind(provision_main_all, provision_main)
provision_main_all_table <- rbind(provision_main_all_table, provision_main_table) %>% 
  select(M_eth18_2011, `n(%, 95%CI)`)

write_csv(provision_main_all, file = "~/provision_main_all.csv")
write_csv(provision_main_all_table, file = "~/provision_main_all_table.csv")

# Provision of baby check by year: main  ----------------------------

# create denominator for each year
denom_main_all_year <- babycheck_mainanalysis %>% 
  group_by(financialyear) %>% 
  count() %>% 
  rename(N = n)

# count number who have received baby check each year
main_all_year <- babycheck_mainanalysis %>% 
  group_by(financialyear, babycheck) %>% 
  count() %>% 
  ungroup() %>% 
  filter(babycheck == 1) %>% 
  select(-babycheck) 

# join and create percentages
main_all_year <- main_all_year %>% 
  left_join(denom_main_all_year, by = "financialyear") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  mutate(M_eth18_2011 = "All")
main_all_year_table <- main_all_year %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, financialyear, `n (%, 95%CI)`) 

# create denominators for each year - by ethnicity
eth_denom_main_all_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
eth_main_all_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear, babycheck) %>% 
  count() %>% 
  ungroup() %>% 
  filter(babycheck == 1) %>% 
  select(-babycheck) 

# join and create percentages
eth_main_all_year <- eth_main_all_year %>% 
  left_join(eth_denom_main_all_year, by = c("M_eth18_2011", "financialyear")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) 
eth_main_all_year_table <- eth_main_all_year %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, financialyear, `n (%, 95%CI)`) 

eth_main_all_year <- rbind(main_all_year, eth_main_all_year) 
eth_main_all_year_table <- rbind(main_all_year_table, eth_main_all_year_table) %>% 
  select(M_eth18_2011,  financialyear,`n (%, 95%CI)`)

write_csv(eth_main_all_year, file="~/eth_main_all_year.csv") 
write_csv(eth_main_all_year_table, file="~/eth_main_all_year_table.csv") 

# relabel ethnicity and set factor level
eth_main_all_year <- eth_main_all_year %>% 
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) %>% 
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "All",
                                    "English, Welsh, Scottish, Northern Irish or British", "Irish", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Any other ethnic group", "Unknown")) %>% 
  rename(`Mother's ethnicity` = M_eth18_2011) 

# dot plots 
eth_main_all_year_dotplot <- ggplot(eth_main_all_year, aes(x = financialyear, y = percent, colour = `Mother's ethnicity`)) + 
  geom_line(aes(group = `Mother's ethnicity`), size = 0.5) +
  theme_bw() +
  geom_point(aes(shape = `Mother's ethnicity`, colour = `Mother's ethnicity`), size = 4) +
  scale_shape_manual(values = c(16, 8, 17, 0, 16, 8, 17, 0, 18, 16, 8, 17,16, 8, 17, 0, 16, 8)) +
  scale_color_manual(values = c("All" = "#564D50",
                                "English, Welsh, Scottish, Northern Irish or British" = "#970055",
                                "Irish" = "#970055",
                                "Any other White background" = "#970055",
                                "Indian" = "#F5B521",
                                "Pakistani" = "#F5B521",
                                "Bangladeshi" = "#F5B521",
                                "Chinese" = "#F5B521",
                                "Any other Asian background" = "#F5B521",
                                "Caribbean" = "#158754",
                                "African" = "#158754",
                                "Any other Black, African or Caribbean background" = "#158754",
                                "White and Black Caribbean" = "#2F84FF",
                                "White and Black African" = "#2F84FF",
                                "White and Asian" = "#2F84FF" ,
                                "Any other Mixed or multiple ethnic background" = "#2F84FF",
                                "Any other ethnic group" = "#A438CD",
                                "Unknown" = "#564D50")) +
  theme(text = element_text(size =20)) +
  labs(y = "Received 6-8 week baby check (%)", x = "Year of birth (financial year)") +
  ylim(50, 100)
eth_main_all_year_dotplot

ggsave(width = 30, height = 8, dpi = 450, "~/eth_main_all_year_dotplot.jpg")

# Provision of baby check by region: main  -------------------------------------------------

# create denominator for each region
denom_main_all_region <- babycheck_mainanalysis %>% 
  group_by(region) %>% 
  count() %>% 
  rename(N = n)

# count number who have received baby check each region
main_all_region <- babycheck_mainanalysis %>% 
  group_by(region, babycheck) %>% 
  count() %>% 
  ungroup() %>% 
  filter(babycheck == 1) %>% 
  select(-babycheck) 

# join and create percentages
main_all_region <- main_all_region %>% 
  left_join(denom_main_all_region, by = "region") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  mutate(M_eth18_2011 = "All")
main_all_region_table <- main_all_region %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, region, `n (%, 95%CI)`) 

# create denominators for each region - by ethnicity
eth_denom_main_all_region <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, region) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each region
eth_main_all_region <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, region, babycheck) %>% 
  count() %>% 
  ungroup() %>% 
  filter(babycheck == 1) %>% 
  select(-babycheck) 

# join and create percentages
eth_main_all_region <- eth_main_all_region %>% 
  left_join(eth_denom_main_all_region, by = c("M_eth18_2011", "region")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  mutate(upperci = replace(upperci, upperci >100, 100))
eth_main_all_region_table <- eth_main_all_region %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, region, `n (%, 95%CI)`) 

eth_main_all_region <- rbind(main_all_region, eth_main_all_region) 
eth_main_all_region_table <- rbind(main_all_region_table, eth_main_all_region_table) %>% 
  select(M_eth18_2011,  region,`n (%, 95%CI)`)

write_csv(eth_main_all_region, file="~/eth_main_all_region.csv") 
write_csv(eth_main_all_region_table, file="~/eth_main_all_region_table.csv") 

# relabel ethnicity and set factor level
eth_main_all_region <- eth_main_all_region %>% 
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) %>% 
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "All",
                                    "English, Welsh, Scottish, Northern Irish or British", "Irish", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Any other ethnic group", "Unknown")) %>% 
  rename(`Mother's ethnicity` = M_eth18_2011) 

# filter out regions were groups are <50
eth_main_all_region <- eth_main_all_region %>% 
  filter(region != "North East") %>% 
  filter(region != "East Midlands") %>% 
  filter(region != "Yorkshire and The Humber")

# dot plots 
pd <- position_dodge(0.7) # move them .07 to the left and right
eth_main_all_region_dotplot <- ggplot(eth_main_all_region, aes(x = region, y = percent, colour = `Mother's ethnicity`)) + 
  geom_errorbar(aes(ymin=lowerci, ymax=upperci), width=.1, position=pd) +
  theme_bw() +
  geom_point(aes(shape = `Mother's ethnicity`, colour = `Mother's ethnicity`), position=pd, size = 4) +
  scale_shape_manual(values = c(16, 8, 17, 0, 16, 8, 17, 0, 18, 16, 8, 17,16, 8, 17, 0, 16, 8)) +
  scale_color_manual(values = c("All" = "#564D50",
                                "English, Welsh, Scottish, Northern Irish or British" = "#970055",
                                "Irish" = "#970055",
                                "Any other White background" = "#970055",
                                "Indian" = "#F5B521",
                                "Pakistani" = "#F5B521",
                                "Bangladeshi" = "#F5B521",
                                "Chinese" = "#F5B521",
                                "Any other Asian background" = "#F5B521",
                                "Caribbean" = "#158754",
                                "African" = "#158754",
                                "Any other Black, African or Caribbean background" = "#158754",
                                "White and Black Caribbean" = "#2F84FF",
                                "White and Black African" = "#2F84FF",
                                "White and Asian" = "#2F84FF" ,
                                "Any other Mixed or multiple ethnic background" = "#2F84FF",
                                "Any other ethnic group" = "#A438CD",
                                "Unknown" = "#564D50")) +
  theme(text = element_text(size =20)) +  labs(y = "Received 6-8 week baby check (%)", x = "Region") +
  ylim(50, 100)
eth_main_all_region_dotplot

ggsave(width = 30, height = 8, dpi = 450, "~/eth_main_all_region_dotplot.jpg")

# Provision of baby check: sens -------------------------------------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

# Provision of baby check by year: sens  ----------------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

# Provision of baby check by region: sens  -------------------------------------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

#--------------------------------------------------------------------------
# UNIVARIABLE MODELLING + EFFECT MODIFICATION                      
#--------------------------------------------------------------------------

load(file="~/babycheck_mainanalysis.Rdata")
load(file="~/babycheck_sensanalysis.Rdata")

# Univariable + regionEM: main ----------------------

# filter out regions were groups are <50
babycheck_mainanalysis <- babycheck_mainanalysis %>% 
  filter(region != "North East") %>% 
  filter(region != "East Midlands") %>% 
  filter(region != "Yorkshire and The Humber")

# run comparator and EM poisson
model_univar_region_mpoisson_main <- glm(babycheck ~ M_eth18_2011 + region, 
                                         data = babycheck_mainanalysis, family = poisson(link = "log"))
model_univar_regionEM_mpoisson_main <- glm(babycheck ~ M_eth18_2011*region, 
                                           data = babycheck_mainanalysis, family = poisson(link = "log"))

# Likelihood ratio test 
lrtest(model_univar_region_mpoisson_main, model_univar_regionEM_mpoisson_main) 

# check p values of interaction term
summary(model_univar_regionEM_mpoisson_main)

# run modified poisson with robust sandwich estimator
univar_mpoisson_regionEM_main_sandwich <- coeftest(model_univar_regionEM_mpoisson_main, vcov = sandwich) 

# extract results and save
univar_mpoisson_regionEM_main_sandwich_clean <- round(cbind(exp(cbind(RR = univar_mpoisson_regionEM_main_sandwich[,1],
                                                                     LCI = univar_mpoisson_regionEM_main_sandwich[,1] + qnorm(0.05/2)*univar_mpoisson_regionEM_main_sandwich[,2],
                                                                     UCI = univar_mpoisson_regionEM_main_sandwich[,1] - qnorm(0.05/2)*univar_mpoisson_regionEM_main_sandwich[,2])),
                                                           P = univar_mpoisson_regionEM_main_sandwich[,4]),2)
univar_mpoisson_regionEM_main_sandwich_clean <- as.data.frame(univar_mpoisson_regionEM_main_sandwich_clean)
extract_mpoisson_results(model_univar_regionEM_mpoisson_main, univar_mpoisson_regionEM_main_sandwich_clean, 'univar_mpoisson_regionEM_main', babycheck_mainanalysis, M_eth18_2011)
write_csv(univar_mpoisson_regionEM_main_table, file="~/modelling/univar_mpoisson_regionEM_main_table.csv") 

# for plotting
write_csv(univar_mpoisson_regionEM_main, file="~/modelling/univar_mpoisson_regionEM_main.csv") 


# Calculate stratum specific effects from unadjusted EM model: main --------------------------------------

ethnicity <- c("Irish", "Any other White background",
               "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
               "Caribbean", "African", "Any other Black, African or Caribbean background",
               "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
               "Any other ethnic group", "Unknown")

## London
univar_mpoisson_regionEM_main_stratum_london <- univar_mpoisson_regionEM_main_sandwich_clean[2:17,] %>% 
  cbind(ethnicity) %>% 
  select(-P) %>% 
  rename(estimate = RR) %>% 
  rename(upper = UCI) %>% 
  rename(lower = LCI) %>% 
  mutate(region = "London")

univar_mpoisson_regionEM_main_stratum_london_table <- univar_mpoisson_regionEM_main_stratum_london %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## East of England

# calculate stratum specific effects - original model
irish_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionEast of England' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionEast of England")
aow_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionEast of England")
indian_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionEast of England")
pakistani_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionEast of England")
bangladeshi_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionEast of England")
chinese_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionEast of England")
aoa_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionEast of England")
caribbean_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionEast of England")
african_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionEast of England")
aob_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionEast of England")
wbc_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionEast of England")
wba_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionEast of England")
wa_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionEast of England")
aom_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionEast of England")
aoe_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionEast of England")
unknown_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionEast of England")

all_stratum_estimates <- rbind(irish_vs_ref_eastofengland, aow_vs_ref_eastofengland,
                               indian_vs_ref_eastofengland, pakistani_vs_ref_eastofengland,
                               bangladeshi_vs_ref_eastofengland, chinese_vs_ref_eastofengland,
                               aoa_vs_ref_eastofengland, caribbean_vs_ref_eastofengland,
                               african_vs_ref_eastofengland, aob_vs_ref_eastofengland,
                               wbc_vs_ref_eastofengland, wba_vs_ref_eastofengland,
                               wa_vs_ref_eastofengland, aom_vs_ref_eastofengland,
                               aoe_vs_ref_eastofengland, unknown_vs_ref_eastofengland)
univar_mpoisson_regionEM_main_stratum_eastofengland <- as.data.frame(t(cbind(data.frame(irish_vs_ref_eastofengland_sandwich), data.frame(aow_vs_ref_eastofengland_sandwich),
                                                                    data.frame(indian_vs_ref_eastofengland_sandwich), data.frame(pakistani_vs_ref_eastofengland_sandwich),
                                                                    data.frame(bangladeshi_vs_ref_eastofengland_sandwich), data.frame(chinese_vs_ref_eastofengland_sandwich),
                                                                    data.frame(aoa_vs_ref_eastofengland_sandwich), data.frame(caribbean_vs_ref_eastofengland_sandwich),
                                                                    data.frame(african_vs_ref_eastofengland_sandwich), data.frame(aob_vs_ref_eastofengland_sandwich),
                                                                    data.frame(wbc_vs_ref_eastofengland_sandwich), data.frame(wba_vs_ref_eastofengland_sandwich),
                                                                    data.frame(wa_vs_ref_eastofengland_sandwich), data.frame(aom_vs_ref_eastofengland_sandwich),
                                                                    data.frame(aoe_vs_ref_eastofengland_sandwich), data.frame(unknown_vs_ref_eastofengland_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "East of England")

univar_mpoisson_regionEM_main_stratum_eastofengland_table <- univar_mpoisson_regionEM_main_stratum_eastofengland %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## North West

# calculate stratum specific effects - original model
irish_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_northwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionNorth West' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionNorth West")
aow_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionNorth West")
indian_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionNorth West")
pakistani_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionNorth West")
bangladeshi_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionNorth West")
chinese_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionNorth West")
aoa_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionNorth West")
caribbean_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionNorth West")
african_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionNorth West")
aob_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionNorth West")
wbc_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionNorth West")
wba_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionNorth West")
wa_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionNorth West")
aom_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionNorth West")
aoe_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionNorth West")
unknown_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionNorth West")

all_stratum_estimates <- rbind(irish_vs_ref_northwest, aow_vs_ref_northwest,
                               indian_vs_ref_northwest, pakistani_vs_ref_northwest,
                               bangladeshi_vs_ref_northwest, chinese_vs_ref_northwest,
                               aoa_vs_ref_northwest, caribbean_vs_ref_northwest,
                               african_vs_ref_northwest, aob_vs_ref_northwest,
                               wbc_vs_ref_northwest, wba_vs_ref_northwest,
                               wa_vs_ref_northwest, aom_vs_ref_northwest,
                               aoe_vs_ref_northwest, unknown_vs_ref_northwest)
univar_mpoisson_regionEM_main_stratum_northwest <- as.data.frame(t(cbind(data.frame(irish_vs_ref_northwest_sandwich), data.frame(aow_vs_ref_northwest_sandwich),
                                                                    data.frame(indian_vs_ref_northwest_sandwich), data.frame(pakistani_vs_ref_northwest_sandwich),
                                                                    data.frame(bangladeshi_vs_ref_northwest_sandwich), data.frame(chinese_vs_ref_northwest_sandwich),
                                                                    data.frame(aoa_vs_ref_northwest_sandwich), data.frame(caribbean_vs_ref_northwest_sandwich),
                                                                    data.frame(african_vs_ref_northwest_sandwich), data.frame(aob_vs_ref_northwest_sandwich),
                                                                    data.frame(wbc_vs_ref_northwest_sandwich), data.frame(wba_vs_ref_northwest_sandwich),
                                                                    data.frame(wa_vs_ref_northwest_sandwich), data.frame(aom_vs_ref_northwest_sandwich),
                                                                    data.frame(aoe_vs_ref_northwest_sandwich), data.frame(unknown_vs_ref_northwest_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "North West")

univar_mpoisson_regionEM_main_stratum_northwest_table <- univar_mpoisson_regionEM_main_stratum_northwest %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## South East

# calculate stratum specific effects - original model
irish_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_southeast <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionSouth East' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionSouth East")
aow_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionSouth East")
indian_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionSouth East")
pakistani_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionSouth East")
bangladeshi_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionSouth East")
chinese_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionSouth East")
aoa_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionSouth East")
caribbean_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionSouth East")
african_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionSouth East")
aob_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionSouth East")
wbc_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionSouth East")
wba_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionSouth East")
wa_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionSouth East")
aom_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth East")
aoe_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionSouth East")
unknown_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionSouth East")

all_stratum_estimates <- rbind(irish_vs_ref_southeast, aow_vs_ref_southeast,
                               indian_vs_ref_southeast, pakistani_vs_ref_southeast,
                               bangladeshi_vs_ref_southeast, chinese_vs_ref_southeast,
                               aoa_vs_ref_southeast, caribbean_vs_ref_southeast,
                               african_vs_ref_southeast, aob_vs_ref_southeast,
                               wbc_vs_ref_southeast, wba_vs_ref_southeast,
                               wa_vs_ref_southeast, aom_vs_ref_southeast,
                               aoe_vs_ref_southeast, unknown_vs_ref_southeast)
univar_mpoisson_regionEM_main_stratum_southeast <- as.data.frame(t(cbind(data.frame(irish_vs_ref_southeast_sandwich), data.frame(aow_vs_ref_southeast_sandwich),
                                                                    data.frame(indian_vs_ref_southeast_sandwich), data.frame(pakistani_vs_ref_southeast_sandwich),
                                                                    data.frame(bangladeshi_vs_ref_southeast_sandwich), data.frame(chinese_vs_ref_southeast_sandwich),
                                                                    data.frame(aoa_vs_ref_southeast_sandwich), data.frame(caribbean_vs_ref_southeast_sandwich),
                                                                    data.frame(african_vs_ref_southeast_sandwich), data.frame(aob_vs_ref_southeast_sandwich),
                                                                    data.frame(wbc_vs_ref_southeast_sandwich), data.frame(wba_vs_ref_southeast_sandwich),
                                                                    data.frame(wa_vs_ref_southeast_sandwich), data.frame(aom_vs_ref_southeast_sandwich),
                                                                    data.frame(aoe_vs_ref_southeast_sandwich), data.frame(unknown_vs_ref_southeast_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "South East")

univar_mpoisson_regionEM_main_stratum_southeast_table <- univar_mpoisson_regionEM_main_stratum_southeast %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## South West

# calculate stratum specific effects - original model
irish_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_southwest <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionSouth West' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionSouth West")
aow_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionSouth West")
indian_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionSouth West")
pakistani_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionSouth West")
bangladeshi_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionSouth West")
chinese_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionSouth West")
aoa_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionSouth West")
caribbean_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionSouth West")
african_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionSouth West")
aob_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionSouth West")
wbc_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionSouth West")
wba_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionSouth West")
wa_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionSouth West")
aom_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth West")
aoe_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionSouth West")
unknown_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionSouth West")

all_stratum_estimates <- rbind(irish_vs_ref_southwest, aow_vs_ref_southwest,
                               indian_vs_ref_southwest, pakistani_vs_ref_southwest,
                               bangladeshi_vs_ref_southwest, chinese_vs_ref_southwest,
                               aoa_vs_ref_southwest, caribbean_vs_ref_southwest,
                               african_vs_ref_southwest, aob_vs_ref_southwest,
                               wbc_vs_ref_southwest, wba_vs_ref_southwest,
                               wa_vs_ref_southwest, aom_vs_ref_southwest,
                               aoe_vs_ref_southwest, unknown_vs_ref_southwest)
univar_mpoisson_regionEM_main_stratum_southwest <- as.data.frame(t(cbind(data.frame(irish_vs_ref_southwest_sandwich), data.frame(aow_vs_ref_southwest_sandwich),
                                                                    data.frame(indian_vs_ref_southwest_sandwich), data.frame(pakistani_vs_ref_southwest_sandwich),
                                                                    data.frame(bangladeshi_vs_ref_southwest_sandwich), data.frame(chinese_vs_ref_southwest_sandwich),
                                                                    data.frame(aoa_vs_ref_southwest_sandwich), data.frame(caribbean_vs_ref_southwest_sandwich),
                                                                    data.frame(african_vs_ref_southwest_sandwich), data.frame(aob_vs_ref_southwest_sandwich),
                                                                    data.frame(wbc_vs_ref_southwest_sandwich), data.frame(wba_vs_ref_southwest_sandwich),
                                                                    data.frame(wa_vs_ref_southwest_sandwich), data.frame(aom_vs_ref_southwest_sandwich),
                                                                    data.frame(aoe_vs_ref_southwest_sandwich), data.frame(unknown_vs_ref_southwest_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "South West")

univar_mpoisson_regionEM_main_stratum_southwest_table <- univar_mpoisson_regionEM_main_stratum_southwest %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## West Midlands

# calculate stratum specific effects - original model
irish_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_univar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionWest Midlands")
aow_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionWest Midlands")
indian_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionWest Midlands")
pakistani_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionWest Midlands")
bangladeshi_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionWest Midlands")
chinese_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionWest Midlands")
aoa_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionWest Midlands")
caribbean_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionWest Midlands")
african_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionWest Midlands")
aob_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionWest Midlands")
wbc_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionWest Midlands")
wba_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionWest Midlands")
wa_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionWest Midlands")
aom_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionWest Midlands")
aoe_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionWest Midlands")
unknown_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionWest Midlands")

all_stratum_estimates <- rbind(irish_vs_ref_westmidlands, aow_vs_ref_westmidlands,
                               indian_vs_ref_westmidlands, pakistani_vs_ref_westmidlands,
                               bangladeshi_vs_ref_westmidlands, chinese_vs_ref_westmidlands,
                               aoa_vs_ref_westmidlands, caribbean_vs_ref_westmidlands,
                               african_vs_ref_westmidlands, aob_vs_ref_westmidlands,
                               wbc_vs_ref_westmidlands, wba_vs_ref_westmidlands,
                               wa_vs_ref_westmidlands, aom_vs_ref_westmidlands,
                               aoe_vs_ref_westmidlands, unknown_vs_ref_westmidlands)
univar_mpoisson_regionEM_main_stratum_westmidlands <- as.data.frame(t(cbind(data.frame(irish_vs_ref_westmidlands_sandwich), data.frame(aow_vs_ref_westmidlands_sandwich),
                                                                  data.frame(indian_vs_ref_westmidlands_sandwich), data.frame(pakistani_vs_ref_westmidlands_sandwich),
                                                                  data.frame(bangladeshi_vs_ref_westmidlands_sandwich), data.frame(chinese_vs_ref_westmidlands_sandwich),
                                                                  data.frame(aoa_vs_ref_westmidlands_sandwich), data.frame(caribbean_vs_ref_westmidlands_sandwich),
                                                                  data.frame(african_vs_ref_westmidlands_sandwich), data.frame(aob_vs_ref_westmidlands_sandwich),
                                                                  data.frame(wbc_vs_ref_westmidlands_sandwich), data.frame(wba_vs_ref_westmidlands_sandwich),
                                                                  data.frame(wa_vs_ref_westmidlands_sandwich), data.frame(aom_vs_ref_westmidlands_sandwich),
                                                                  data.frame(aoe_vs_ref_westmidlands_sandwich), data.frame(unknown_vs_ref_westmidlands_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "West Midlands")

univar_mpoisson_regionEM_main_stratum_westmidlands_table <- univar_mpoisson_regionEM_main_stratum_westmidlands %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

# combine regions
univar_mpoisson_regionEM_main_strata <- rbind(univar_mpoisson_regionEM_main_stratum_london,
                                              univar_mpoisson_regionEM_main_stratum_eastofengland,
                                              univar_mpoisson_regionEM_main_stratum_northwest,
                                              univar_mpoisson_regionEM_main_stratum_southeast,
                                              univar_mpoisson_regionEM_main_stratum_southwest,
                                            univar_mpoisson_regionEM_main_stratum_westmidlands)

univar_mpoisson_regionEM_main_strata_table <-rbind(univar_mpoisson_regionEM_main_stratum_london_table,
                                                   univar_mpoisson_regionEM_main_stratum_eastofengland_table,
                                                   univar_mpoisson_regionEM_main_stratum_northwest_table,
                                                   univar_mpoisson_regionEM_main_stratum_southeast_table,
                                                   univar_mpoisson_regionEM_main_stratum_southwest_table,
                                                 univar_mpoisson_regionEM_main_stratum_westmidlands_table)

# save
write_csv(univar_mpoisson_regionEM_main_strata, file="~/modelling/univar_mpoisson_regionEM_main_strata.csv") 
write_csv(univar_mpoisson_regionEM_main_strata_table, file="~/modelling/univar_mpoisson_regionEM_main_strata_table.csv") 

# remove model
rm(model_univar_region_mpoisson_main, model_univar_regionEM_mpoisson_main)


# Univariable + regionEM: sens ----------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

# Calculate stratum specific effects from unadjusted EM model: sens --------------------------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed


#--------------------------------------------------------------------------
# MULTIVARIABLE MODELLING + EFFECT MODIFICATION                    
#--------------------------------------------------------------------------

load(file="~/babycheck_mainanalysis.Rdata")
load(file="~/babycheck_sensanalysis.Rdata")

# Set reference groups ----------------------------------------------------

babycheck_mainanalysis$matage <- relevel(babycheck_mainanalysis$matage, ref = "30-34") # as per Li et al maternal checks paper
babycheck_mainanalysis$gest <- relevel(babycheck_mainanalysis$gest, ref = "term (37-42 weeks)")
babycheck_sensanalysis$matage <- relevel(babycheck_sensanalysis$matage, ref = "30-34")
babycheck_sensanalysis$gest <- relevel(babycheck_sensanalysis$gest, ref = "term (37-42 weeks)")

# Adjusted for SD + mat/birth mediators + regionEM: main  ----------------------

babycheck_mainanalysis <- babycheck_mainanalysis %>% 
  filter(region != "North East") %>% 
  filter(region != "East Midlands") %>% 
  filter(region != "Yorkshire and The Humber")

# run comparator and EM poisson
model_multivar_region_mpoisson_main <- glm(babycheck ~ M_eth18_2011 + region +
                                           imd + ruc +
                                           preterm + modebirth + matage + first_time_mum + 
                                           hospitalised, 
                                         data = babycheck_mainanalysis, family = poisson(link = "log"))
model_multivar_regionEM_mpoisson_main <- glm(babycheck ~ M_eth18_2011*region +
                                             imd + ruc +
                                             preterm + modebirth + matage + first_time_mum + 
                                             hospitalised, 
                                           data = babycheck_mainanalysis, family = poisson(link = "log"))

# Likelihood ratio test 
lrtest(model_multivar_region_mpoisson_main, model_multivar_regionEM_mpoisson_main) 

# check p values of interaction term
summary(model_multivar_regionEM_mpoisson_main)

# run modified poisson with robust sandwich estimator
multivar_mpoisson_regionEM_main_sandwich <- coeftest(model_multivar_regionEM_mpoisson_main, vcov = sandwich) 

# extract results and save
multivar_mpoisson_regionEM_main_sandwich_clean <- round(cbind(exp(cbind(RR = multivar_mpoisson_regionEM_main_sandwich[,1],
                                                                      LCI = multivar_mpoisson_regionEM_main_sandwich[,1] + qnorm(0.05/2)*multivar_mpoisson_regionEM_main_sandwich[,2],
                                                                      UCI = multivar_mpoisson_regionEM_main_sandwich[,1] - qnorm(0.05/2)*multivar_mpoisson_regionEM_main_sandwich[,2])),
                                                            P = multivar_mpoisson_regionEM_main_sandwich[,4]),2)
multivar_mpoisson_regionEM_main_sandwich_clean <- as.data.frame(multivar_mpoisson_regionEM_main_sandwich_clean)
extract_mpoisson_results(model_multivar_regionEM_mpoisson_main, multivar_mpoisson_regionEM_main_sandwich_clean, 'multivar_mpoisson_regionEM_main', babycheck_mainanalysis, M_eth18_2011)
write_csv(multivar_mpoisson_regionEM_main_table, file="~/modelling/multivar_mpoisson_regionEM_main_table.csv") 

# for plotting
write_csv(multivar_mpoisson_regionEM_main, file="~/modelling/multivar_mpoisson_regionEM_main.csv") 

# Calculate stratum specific effects from adjusted EM model: main --------------------------------------

ethnicity <- c("Irish", "Any other White background",
               "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
               "Caribbean", "African", "Any other Black, African or Caribbean background",
               "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
               "Any other ethnic group", "Unknown")

## London
multivar_mpoisson_regionEM_main_stratum_london <- multivar_mpoisson_regionEM_main_sandwich_clean[2:17,] %>% 
  cbind(ethnicity) %>% 
  select(-P) %>% 
  rename(estimate = RR) %>% 
  rename(upper = UCI) %>% 
  rename(lower = LCI) %>% 
  mutate(region = "London")

multivar_mpoisson_regionEM_main_stratum_london_table <- multivar_mpoisson_regionEM_main_stratum_london %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## East of England

# calculate stratum specific effects - original model
irish_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionEast of England' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_eastofengland <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionEast of England' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionEast of England")
aow_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionEast of England")
indian_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionEast of England")
pakistani_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionEast of England")
bangladeshi_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionEast of England")
chinese_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionEast of England")
aoa_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionEast of England")
caribbean_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionEast of England")
african_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionEast of England")
aob_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionEast of England")
wbc_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionEast of England")
wba_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionEast of England")
wa_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionEast of England")
aom_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionEast of England")
aoe_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionEast of England")
unknown_vs_ref_eastofengland_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionEast of England")

all_stratum_estimates <- rbind(irish_vs_ref_eastofengland, aow_vs_ref_eastofengland,
                               indian_vs_ref_eastofengland, pakistani_vs_ref_eastofengland,
                               bangladeshi_vs_ref_eastofengland, chinese_vs_ref_eastofengland,
                               aoa_vs_ref_eastofengland, caribbean_vs_ref_eastofengland,
                               african_vs_ref_eastofengland, aob_vs_ref_eastofengland,
                               wbc_vs_ref_eastofengland, wba_vs_ref_eastofengland,
                               wa_vs_ref_eastofengland, aom_vs_ref_eastofengland,
                               aoe_vs_ref_eastofengland, unknown_vs_ref_eastofengland)
multivar_mpoisson_regionEM_main_stratum_eastofengland <- as.data.frame(t(cbind(data.frame(irish_vs_ref_eastofengland_sandwich), data.frame(aow_vs_ref_eastofengland_sandwich),
                                                                             data.frame(indian_vs_ref_eastofengland_sandwich), data.frame(pakistani_vs_ref_eastofengland_sandwich),
                                                                             data.frame(bangladeshi_vs_ref_eastofengland_sandwich), data.frame(chinese_vs_ref_eastofengland_sandwich),
                                                                             data.frame(aoa_vs_ref_eastofengland_sandwich), data.frame(caribbean_vs_ref_eastofengland_sandwich),
                                                                             data.frame(african_vs_ref_eastofengland_sandwich), data.frame(aob_vs_ref_eastofengland_sandwich),
                                                                             data.frame(wbc_vs_ref_eastofengland_sandwich), data.frame(wba_vs_ref_eastofengland_sandwich),
                                                                             data.frame(wa_vs_ref_eastofengland_sandwich), data.frame(aom_vs_ref_eastofengland_sandwich),
                                                                             data.frame(aoe_vs_ref_eastofengland_sandwich), data.frame(unknown_vs_ref_eastofengland_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "East of England")

multivar_mpoisson_regionEM_main_stratum_eastofengland_table <- multivar_mpoisson_regionEM_main_stratum_eastofengland %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## North West

# calculate stratum specific effects - original model
irish_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionNorth West' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_northwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionNorth West' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionNorth West")
aow_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionNorth West")
indian_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionNorth West")
pakistani_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionNorth West")
bangladeshi_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionNorth West")
chinese_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionNorth West")
aoa_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionNorth West")
caribbean_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionNorth West")
african_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionNorth West")
aob_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionNorth West")
wbc_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionNorth West")
wba_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionNorth West")
wa_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionNorth West")
aom_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionNorth West")
aoe_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionNorth West")
unknown_vs_ref_northwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionNorth West")

all_stratum_estimates <- rbind(irish_vs_ref_northwest, aow_vs_ref_northwest,
                               indian_vs_ref_northwest, pakistani_vs_ref_northwest,
                               bangladeshi_vs_ref_northwest, chinese_vs_ref_northwest,
                               aoa_vs_ref_northwest, caribbean_vs_ref_northwest,
                               african_vs_ref_northwest, aob_vs_ref_northwest,
                               wbc_vs_ref_northwest, wba_vs_ref_northwest,
                               wa_vs_ref_northwest, aom_vs_ref_northwest,
                               aoe_vs_ref_northwest, unknown_vs_ref_northwest)
multivar_mpoisson_regionEM_main_stratum_northwest <- as.data.frame(t(cbind(data.frame(irish_vs_ref_northwest_sandwich), data.frame(aow_vs_ref_northwest_sandwich),
                                                                         data.frame(indian_vs_ref_northwest_sandwich), data.frame(pakistani_vs_ref_northwest_sandwich),
                                                                         data.frame(bangladeshi_vs_ref_northwest_sandwich), data.frame(chinese_vs_ref_northwest_sandwich),
                                                                         data.frame(aoa_vs_ref_northwest_sandwich), data.frame(caribbean_vs_ref_northwest_sandwich),
                                                                         data.frame(african_vs_ref_northwest_sandwich), data.frame(aob_vs_ref_northwest_sandwich),
                                                                         data.frame(wbc_vs_ref_northwest_sandwich), data.frame(wba_vs_ref_northwest_sandwich),
                                                                         data.frame(wa_vs_ref_northwest_sandwich), data.frame(aom_vs_ref_northwest_sandwich),
                                                                         data.frame(aoe_vs_ref_northwest_sandwich), data.frame(unknown_vs_ref_northwest_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "North West")

multivar_mpoisson_regionEM_main_stratum_northwest_table <- multivar_mpoisson_regionEM_main_stratum_northwest %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## South East

# calculate stratum specific effects - original model
irish_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionSouth East' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_southeast <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionSouth East' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionSouth East")
aow_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionSouth East")
indian_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionSouth East")
pakistani_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionSouth East")
bangladeshi_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionSouth East")
chinese_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionSouth East")
aoa_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionSouth East")
caribbean_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionSouth East")
african_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionSouth East")
aob_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionSouth East")
wbc_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionSouth East")
wba_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionSouth East")
wa_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionSouth East")
aom_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth East")
aoe_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionSouth East")
unknown_vs_ref_southeast_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionSouth East")

all_stratum_estimates <- rbind(irish_vs_ref_southeast, aow_vs_ref_southeast,
                               indian_vs_ref_southeast, pakistani_vs_ref_southeast,
                               bangladeshi_vs_ref_southeast, chinese_vs_ref_southeast,
                               aoa_vs_ref_southeast, caribbean_vs_ref_southeast,
                               african_vs_ref_southeast, aob_vs_ref_southeast,
                               wbc_vs_ref_southeast, wba_vs_ref_southeast,
                               wa_vs_ref_southeast, aom_vs_ref_southeast,
                               aoe_vs_ref_southeast, unknown_vs_ref_southeast)
multivar_mpoisson_regionEM_main_stratum_southeast <- as.data.frame(t(cbind(data.frame(irish_vs_ref_southeast_sandwich), data.frame(aow_vs_ref_southeast_sandwich),
                                                                         data.frame(indian_vs_ref_southeast_sandwich), data.frame(pakistani_vs_ref_southeast_sandwich),
                                                                         data.frame(bangladeshi_vs_ref_southeast_sandwich), data.frame(chinese_vs_ref_southeast_sandwich),
                                                                         data.frame(aoa_vs_ref_southeast_sandwich), data.frame(caribbean_vs_ref_southeast_sandwich),
                                                                         data.frame(african_vs_ref_southeast_sandwich), data.frame(aob_vs_ref_southeast_sandwich),
                                                                         data.frame(wbc_vs_ref_southeast_sandwich), data.frame(wba_vs_ref_southeast_sandwich),
                                                                         data.frame(wa_vs_ref_southeast_sandwich), data.frame(aom_vs_ref_southeast_sandwich),
                                                                         data.frame(aoe_vs_ref_southeast_sandwich), data.frame(unknown_vs_ref_southeast_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "South East")

multivar_mpoisson_regionEM_main_stratum_southeast_table <- multivar_mpoisson_regionEM_main_stratum_southeast %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## South West

# calculate stratum specific effects - original model
irish_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionSouth West' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_southwest <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionSouth West' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionSouth West")
aow_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionSouth West")
indian_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionSouth West")
pakistani_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionSouth West")
bangladeshi_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionSouth West")
chinese_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionSouth West")
aoa_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionSouth West")
caribbean_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionSouth West")
african_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionSouth West")
aob_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionSouth West")
wbc_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionSouth West")
wba_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionSouth West")
wa_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionSouth West")
aom_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionSouth West")
aoe_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionSouth West")
unknown_vs_ref_southwest_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionSouth West")

all_stratum_estimates <- rbind(irish_vs_ref_southwest, aow_vs_ref_southwest,
                               indian_vs_ref_southwest, pakistani_vs_ref_southwest,
                               bangladeshi_vs_ref_southwest, chinese_vs_ref_southwest,
                               aoa_vs_ref_southwest, caribbean_vs_ref_southwest,
                               african_vs_ref_southwest, aob_vs_ref_southwest,
                               wbc_vs_ref_southwest, wba_vs_ref_southwest,
                               wa_vs_ref_southwest, aom_vs_ref_southwest,
                               aoe_vs_ref_southwest, unknown_vs_ref_southwest)
multivar_mpoisson_regionEM_main_stratum_southwest <- as.data.frame(t(cbind(data.frame(irish_vs_ref_southwest_sandwich), data.frame(aow_vs_ref_southwest_sandwich),
                                                                         data.frame(indian_vs_ref_southwest_sandwich), data.frame(pakistani_vs_ref_southwest_sandwich),
                                                                         data.frame(bangladeshi_vs_ref_southwest_sandwich), data.frame(chinese_vs_ref_southwest_sandwich),
                                                                         data.frame(aoa_vs_ref_southwest_sandwich), data.frame(caribbean_vs_ref_southwest_sandwich),
                                                                         data.frame(african_vs_ref_southwest_sandwich), data.frame(aob_vs_ref_southwest_sandwich),
                                                                         data.frame(wbc_vs_ref_southwest_sandwich), data.frame(wba_vs_ref_southwest_sandwich),
                                                                         data.frame(wa_vs_ref_southwest_sandwich), data.frame(aom_vs_ref_southwest_sandwich),
                                                                         data.frame(aoe_vs_ref_southwest_sandwich), data.frame(unknown_vs_ref_southwest_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "South West")

multivar_mpoisson_regionEM_main_stratum_southwest_table <- multivar_mpoisson_regionEM_main_stratum_southwest %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

## West Midlands

# calculate stratum specific effects - original model
irish_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011African' = 1, 'M_eth18_2011African:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_westmidlands <- as.data.frame(exp(estimable(model_multivar_regionEM_mpoisson_main, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:regionWest Midlands' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:regionWest Midlands")
aow_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:regionWest Midlands")
indian_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:regionWest Midlands")
pakistani_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:regionWest Midlands")
bangladeshi_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:regionWest Midlands")
chinese_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:regionWest Midlands")
aoa_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:regionWest Midlands")
caribbean_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:regionWest Midlands")
african_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011African", "M_eth18_2011African:regionWest Midlands")
aob_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:regionWest Midlands")
wbc_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:regionWest Midlands")
wba_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:regionWest Midlands")
wa_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:regionWest Midlands")
aom_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:regionWest Midlands")
aoe_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:regionWest Midlands")
unknown_vs_ref_westmidlands_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_regionEM_main_sandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:regionWest Midlands")

all_stratum_estimates <- rbind(irish_vs_ref_westmidlands, aow_vs_ref_westmidlands,
                               indian_vs_ref_westmidlands, pakistani_vs_ref_westmidlands,
                               bangladeshi_vs_ref_westmidlands, chinese_vs_ref_westmidlands,
                               aoa_vs_ref_westmidlands, caribbean_vs_ref_westmidlands,
                               african_vs_ref_westmidlands, aob_vs_ref_westmidlands,
                               wbc_vs_ref_westmidlands, wba_vs_ref_westmidlands,
                               wa_vs_ref_westmidlands, aom_vs_ref_westmidlands,
                               aoe_vs_ref_westmidlands, unknown_vs_ref_westmidlands)
multivar_mpoisson_regionEM_main_stratum_westmidlands <- as.data.frame(t(cbind(data.frame(irish_vs_ref_westmidlands_sandwich), data.frame(aow_vs_ref_westmidlands_sandwich),
                                                                            data.frame(indian_vs_ref_westmidlands_sandwich), data.frame(pakistani_vs_ref_westmidlands_sandwich),
                                                                            data.frame(bangladeshi_vs_ref_westmidlands_sandwich), data.frame(chinese_vs_ref_westmidlands_sandwich),
                                                                            data.frame(aoa_vs_ref_westmidlands_sandwich), data.frame(caribbean_vs_ref_westmidlands_sandwich),
                                                                            data.frame(african_vs_ref_westmidlands_sandwich), data.frame(aob_vs_ref_westmidlands_sandwich),
                                                                            data.frame(wbc_vs_ref_westmidlands_sandwich), data.frame(wba_vs_ref_westmidlands_sandwich),
                                                                            data.frame(wa_vs_ref_westmidlands_sandwich), data.frame(aom_vs_ref_westmidlands_sandwich),
                                                                            data.frame(aoe_vs_ref_westmidlands_sandwich), data.frame(unknown_vs_ref_westmidlands_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(region = "West Midlands")

multivar_mpoisson_regionEM_main_stratum_westmidlands_table <- multivar_mpoisson_regionEM_main_stratum_westmidlands %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, region)

# combine regions
multivar_mpoisson_regionEM_main_strata <- rbind(multivar_mpoisson_regionEM_main_stratum_london,
                                              multivar_mpoisson_regionEM_main_stratum_eastofengland,
                                              multivar_mpoisson_regionEM_main_stratum_northwest,
                                              multivar_mpoisson_regionEM_main_stratum_southeast,
                                              multivar_mpoisson_regionEM_main_stratum_southwest,
                                              multivar_mpoisson_regionEM_main_stratum_westmidlands)

multivar_mpoisson_regionEM_main_strata_table <-rbind(multivar_mpoisson_regionEM_main_stratum_london_table,
                                                   multivar_mpoisson_regionEM_main_stratum_eastofengland_table,
                                                   multivar_mpoisson_regionEM_main_stratum_northwest_table,
                                                   multivar_mpoisson_regionEM_main_stratum_southeast_table,
                                                   multivar_mpoisson_regionEM_main_stratum_southwest_table,
                                                   multivar_mpoisson_regionEM_main_stratum_westmidlands_table)

# save
write_csv(multivar_mpoisson_regionEM_main_strata, file="~/modelling/multivar_mpoisson_regionEM_main_strata.csv") 
write_csv(multivar_mpoisson_regionEM_main_strata_table, file="~/modelling/multivar_mpoisson_regionEM_main_strata_table.csv") 

# remove model
rm(model_multivar_region_mpoisson_main, model_multivar_regionEM_mpoisson_main)

# Adjusted for SD + mat/birth mediators + regionEM: sens  ----------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

# Calculate stratum specific effects from adjusted EM model: sens --------------------------------------

# Repeat for sensitivity analysis using babycheck_sensanalysis file if needed

#--------------------------------------------------------------------------
# FORESTPLOTS                     
#--------------------------------------------------------------------------

# Unadjusted: main---------------------------------------------------

# load 
unadj_babycheck_strata_main <- read_csv(file="~/modelling/univar_mpoisson_regionEM_main_strata.csv") 

# plot pc
dev.new()
png(file = "~/modelling/fp_unadj_babycheck_strata_main.png", width = 1500, height = 2000)
fp_unadj_babycheck_strata <- unadj_babycheck_strata_main %>% 
  group_by(region) %>%
  forestplot(label = ethnicity,
             mean = estimate,
             lower = lower,
             upper = upper,
             fn.ci_norm = c(fpDrawNormalCI,  fpDrawCircleCI, fpDrawPointCI,fpDrawDiamondCI, fpDrawNormalCI, fpDrawCircleCI),
             boxsize = .125, # We set the box size to better visualize the type
             line.margin = .1, # We need to add this to avoid crowding
             xlog = TRUE,
             clip = c(0.1, 1.3), 
             xticks = log10(c(0.44, 0.60, 0.79, 1, 1.24, 1.51)),
             zero = 1,
             xlab = "RR (95%CI)",
             mar = unit(rep(10, times = 4), "mm"),
             col = fpColors(box = c("black", "grey", "black", "black", "grey", "black"),
                            lines = c("black", "grey", "black", "black", "grey", "black"),
                            zero = "black"),
             txt_gp = fpTxtGp(ticks=gpar(cex=2.5), xlab=gpar(cex=2.5), cex = 2.5, summary = gpar(fontface = 'bold'),
                              label = list(gpar(fontface = 'plain', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5)))) %>% 
  fp_set_zebra_style("#F8F8F8") %>% 
  fp_add_lines(h_3 = gpar(lwd = 1.5, col = "#97928C"), 
               h_8 = gpar(lwd = 1.5, col = "#97928C"),
               h_11 = gpar(lwd = 1.5, col = "#97928C"),
               h_15 = gpar(lwd = 1.5, col = "#97928C"),
               h_16 = gpar(lwd = 1.5, col = "#97928C")) %>% 
  fp_decorate_graph(grid = structure(c(0.7, 0.801, 0.903, 1.098), 
                                     gp = gpar(lty = 2, col = "#D3D3D3"))) 
print(fp_unadj_babycheck_strata)
dev.off()

# Unadjusted: sens ---------------------------------------------------

# Repeat for univar_mpoisson_regionEM_sens_strata as necessary

# Adjusted: main---------------------------------------------------

# load 
adj_babycheck_strata_main <- read_csv(file="~/modelling/multivar_mpoisson_regionEM_main_strata.csv") 

# plot pc
dev.new()
png(file = "~/modelling/fp_adj_babycheck_strata_main.png", width = 1500, height = 2000)
fp_adj_babycheck_strata <- adj_babycheck_strata_main %>% 
  group_by(region) %>%
  forestplot(label = ethnicity,
             mean = estimate,
             lower = lower,
             upper = upper,
             fn.ci_norm = c(fpDrawNormalCI,  fpDrawCircleCI, fpDrawPointCI,fpDrawDiamondCI, fpDrawNormalCI, fpDrawCircleCI),
             boxsize = .125, # We set the box size to better visualize the type
             line.margin = .1, # We need to add this to avoid crowding
             xlog = TRUE,
             clip = c(0.1, 1.3),
             xticks = log10(c(0.44, 0.60, 0.79, 1, 1.24, 1.51)), 
             zero = 1,
             xlab = "RR (95%CI)",
             mar = unit(rep(10, times = 4), "mm"),
             col = fpColors(box = c("black", "grey", "black", "black", "grey", "black"),
                            lines = c("black", "grey", "black", "black", "grey", "black"),
                            zero = "black"),
             txt_gp = fpTxtGp(ticks=gpar(cex=2.5), xlab=gpar(cex=2.5), cex = 2.5, summary = gpar(fontface = 'bold'),
                              label = list(gpar(fontface = 'plain', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5)))) %>% 
  fp_set_zebra_style("#F8F8F8") %>% 
  fp_add_lines(h_3 = gpar(lwd = 1.5, col = "#97928C"), 
               h_8 = gpar(lwd = 1.5, col = "#97928C"),
               h_11 = gpar(lwd = 1.5, col = "#97928C"),
               h_15 = gpar(lwd = 1.5, col = "#97928C"),
               h_16 = gpar(lwd = 1.5, col = "#97928C")) %>% 
  fp_decorate_graph(grid = structure(c(0.7, 0.801, 0.903, 1.098), 
                                     gp = gpar(lty = 2, col = "#D3D3D3"))) 
print(fp_adj_babycheck_strata)
dev.off()

# Adjusted: sens ---------------------------------------------------

# Repeat for multivar_mpoisson_regionEM_sens_strata as necessary


#--------------------------------------------------------------------------
# DESCRIPTIVE ANALYSIS: SERVICE DELIVERY                      
#--------------------------------------------------------------------------

load(file="~/babycheck_mainanalysis.Rdata")
load(file="~/babycheck_sensanalysis.Rdata")

# Provision of maternal check: main -------------------------------------------------

# main analysis, all checks
eth_denom_main_all <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_maternalcheck_main_all <- babycheck_mainanalysis %>% 
  group_by(maternalcheck_4to12reg, M_eth18_2011) %>% 
  count() %>% 
  ungroup()
provision_maternalcheck_main_all <- provision_maternalcheck_main_all %>% 
  left_join(eth_denom_main_all, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(maternalcheck_4to12reg == 1)  
provision_maternalcheck_main_all_table <- provision_maternalcheck_main_all %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, maternalcheck_4to12reg)) 

write_csv(provision_maternalcheck_main_all, file = "~/provision_maternalcheck_main_all.csv")
write_csv(provision_maternalcheck_main_all_table, file = "~/provision_maternalcheck_main_all_table.csv")

# main analysis, all checks, by year
eth_denom_main_all_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_maternalcheck_main_all_year <- babycheck_mainanalysis %>% 
  group_by(maternalcheck_4to12reg, M_eth18_2011, financialyear) %>% 
  count() %>% 
  ungroup() 
provision_maternalcheck_main_all_year <- provision_maternalcheck_main_all_year %>% 
  left_join(eth_denom_main_all_year, by = c("M_eth18_2011", "financialyear")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(maternalcheck_4to12reg == 1) %>% 
  arrange(financialyear)
provision_maternalcheck_main_all_year_table <- provision_maternalcheck_main_all_year %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, maternalcheck_4to12reg)) 

write_csv(provision_maternalcheck_main_all_year, file = "~/provision_maternalcheck_main_all_year.csv")
write_csv(provision_maternalcheck_main_all_year_table, file = "~/provision_maternalcheck_main_all_year_table.csv")


# Provision of maternal check: sens -------------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Timing of baby check + maternal check: main  -----------------------------------------------

# (1 = maternal check before baby check, 2 = both on the same date, 3 = maternal check after baby check, 4 = no maternal check)

# only include those who had a baby check
babycheck_mainanalysis_yescheck <- babycheck_mainanalysis %>% 
  filter(babycheck == 1)

# main analysis, all checks
eth_denom_main_babycheck_maternalcheck <- babycheck_mainanalysis_yescheck %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_babycheck_maternalcheck <- babycheck_mainanalysis_yescheck %>% 
  group_by(M_eth18_2011, babycheck_maternalcheck_4to12reg) %>% 
  count() %>% 
  ungroup()
provision_main_babycheck_maternalcheck <- provision_main_babycheck_maternalcheck %>% 
  left_join(eth_denom_main_babycheck_maternalcheck, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) 
provision_main_babycheck_maternalcheck_table <- provision_main_babycheck_maternalcheck %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci)) %>% 
  arrange(babycheck_maternalcheck_4to12reg)

write_csv(provision_main_babycheck_maternalcheck, file = "~/provision_main_babycheck_maternalcheck.csv")
write_csv(provision_main_babycheck_maternalcheck_table, file = "~/provision_main_babycheck_maternalcheck_table.csv")

# Timing of baby check + maternal check: sens  -----------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Baby check if have maternal check: main  -----------------------------------------------

# main analysis, all checks
eth_denom_main_maternalcheck <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, maternalcheck_4to12reg) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_maternalcheck <- babycheck_mainanalysis %>% 
  group_by(babycheck, M_eth18_2011, maternalcheck_4to12reg) %>% 
  count() %>% 
  ungroup()
provision_main_maternalcheck <- provision_main_maternalcheck %>% 
  left_join(eth_denom_main_maternalcheck, by = c("M_eth18_2011", "maternalcheck_4to12reg")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(babycheck == 1)
provision_main_maternalcheck_table <- provision_main_maternalcheck %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci)) %>% 
  arrange(maternalcheck_4to12reg)

write_csv(provision_main_maternalcheck, file = "~/provision_main_maternalcheck.csv")
write_csv(provision_main_maternalcheck_table, file = "~/provision_main_maternalcheck_table.csv")

# main analysis, all checks, by year
eth_denom_main_maternalcheck_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear, maternalcheck_4to12reg) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_maternalcheck_year <- babycheck_mainanalysis %>% 
  group_by(babycheck,M_eth18_2011, maternalcheck_4to12reg, financialyear) %>% 
  count() %>% 
  ungroup() 
provision_main_maternalcheck_year <- provision_main_maternalcheck_year %>% 
  left_join(eth_denom_main_maternalcheck_year, by = c("M_eth18_2011", "financialyear","maternalcheck_4to12reg")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  arrange(financialyear) %>% 
  filter(babycheck == 1)
provision_main_maternalcheck_year_table <- provision_main_maternalcheck_year %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, babycheck)) 

write_csv(provision_main_maternalcheck_year, file = "~/provision_main_maternalcheck_year.csv")
write_csv(provision_main_maternalcheck_year_table, file = "~/provision_main_maternalcheck_year_table.csv")

# Baby check if have maternal check: sens -----------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Infant vax: main -------------------------------------------------

# main analysis, all checks
eth_denom_main_all <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_infantvax_main_all <- babycheck_mainanalysis %>% 
  group_by(infantvax, M_eth18_2011) %>% 
  count() %>% 
  ungroup()
provision_infantvax_main_all <- provision_infantvax_main_all %>% 
  left_join(eth_denom_main_all, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(infantvax == 1)  
provision_infantvax_main_all_table <- provision_infantvax_main_all %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, infantvax)) 

write_csv(provision_infantvax_main_all, file = "~/provision_infantvax_main_all.csv")
write_csv(provision_infantvax_main_all_table, file = "~/provision_infantvax_main_all_table.csv")

# main analysis, all checks, by year
eth_denom_main_all_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_infantvax_main_all_year <- babycheck_mainanalysis %>% 
  group_by(infantvax, M_eth18_2011, financialyear) %>% 
  count() %>% 
  ungroup() 
provision_infantvax_main_all_year <- provision_infantvax_main_all_year %>% 
  left_join(eth_denom_main_all_year, by = c("M_eth18_2011", "financialyear")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(infantvax == 1) %>% 
  arrange(financialyear)
provision_infantvax_main_all_year_table <- provision_infantvax_main_all_year %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, infantvax)) 

write_csv(provision_infantvax_main_all_year, file = "~/provision_infantvax_main_all_year.csv")
write_csv(provision_infantvax_main_all_year_table, file = "~/provision_infantvax_main_all_year_table.csv")



# Infant vax: sens -------------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Timing of baby check + vaccination: main -------------------------------------------------------------

# only include those who had a baby check
babycheck_mainanalysis_yescheck <- babycheck_mainanalysis %>% 
  filter(babycheck == 1)

# main analysis, all checks
eth_denom_main_babycheck_vax <- babycheck_mainanalysis_yescheck %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_babycheck_vax <- babycheck_mainanalysis_yescheck %>% 
  group_by(babycheck, M_eth18_2011, babycheck_vax) %>% 
  count() %>% 
  ungroup()
provision_main_babycheck_vax <- provision_main_babycheck_vax %>% 
  left_join(eth_denom_main_babycheck_vax, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(babycheck_vax == 1)
provision_main_babycheck_vax_table <- provision_main_babycheck_vax %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, babycheck)) %>% 
  arrange(babycheck_vax)

write_csv(provision_main_babycheck_vax, file = "~/provision_main_babycheck_vax.csv")
write_csv(provision_main_babycheck_vax_table, file = "~/provision_main_babycheck_vax_table.csv")

# Timing of baby check + vaccination: sens -------------------------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Baby check if have infant vax: main -----------------------------------------------

# main analysis, all checks
eth_denom_main_infantvax <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, infantvax) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_infantvax <- babycheck_mainanalysis %>% 
  group_by(babycheck, M_eth18_2011, infantvax) %>% 
  count() %>% 
  ungroup()
provision_main_infantvax <- provision_main_infantvax %>% 
  left_join(eth_denom_main_infantvax, by = c("M_eth18_2011", "infantvax")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(babycheck == 1)
provision_main_infantvax_table <- provision_main_infantvax %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci)) %>% 
  arrange(infantvax)

write_csv(provision_main_infantvax, file = "~/provision_main_infantvax.csv")
write_csv(provision_main_infantvax_table, file = "~/provision_main_infantvax_table.csv")

# main analysis, all checks, by year
eth_denom_main_infantvax_year <- babycheck_mainanalysis %>% 
  group_by(M_eth18_2011, financialyear, infantvax) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
provision_main_infantvax_year <- babycheck_mainanalysis %>% 
  group_by(babycheck,M_eth18_2011, infantvax, financialyear) %>% 
  count() %>% 
  ungroup() 
provision_main_infantvax_year <- provision_main_infantvax_year %>% 
  left_join(eth_denom_main_infantvax_year, by = c("M_eth18_2011", "financialyear","infantvax")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  arrange(financialyear) %>% 
  filter(babycheck == 1)
provision_main_infantvax_year_table <- provision_main_infantvax_year %>% 
  mutate(`n(%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  dplyr::select(-c(n, N, percent, prop, upperci,lowerci, babycheck)) 

write_csv(provision_main_infantvax_year, file = "~/provision_main_infantvax_year.csv")
write_csv(provision_main_infantvax_year_table, file = "~/provision_main_infantvax_year_table.csv")

# Baby check if have infant vax: sens  -----------------------------------------------

# Repeat for babycheck_sensanalysis as necessary 

# Baby check + maternal check + infant vax: main --------------------------------

# create new variable for who had all 3 services on the same day vs 2 services vs just baby check
# 1 = just baby check on that day, 2 = baby check + maternal check on the same day, 3 = baby check + vax on the same day, 4 = all 3 on the same day 
# only include those who had a baby check
babycheck_mainanalysis_allthree <- babycheck_mainanalysis %>% 
  filter(babycheck == 1)

babycheck_mainanalysis_allthree <- babycheck_mainanalysis_allthree %>% 
  mutate(allthree = "Baby check") %>% 
  mutate(allthree = replace(allthree, babycheck_vax == 0 & babycheck_maternalcheck_4to12reg == 2, "Baby check & maternal check")) %>%  
  mutate(allthree = replace(allthree, babycheck_vax == 1 & babycheck_maternalcheck_4to12reg != 2, "Baby check & infant vaccination")) %>% 
  mutate(allthree = replace(allthree, babycheck_vax == 1 & babycheck_maternalcheck_4to12reg == 2, "All 3 services")) %>% 
  mutate(allthree = as.factor(allthree)) %>% 
  mutate(allthree = fct_relevel(allthree, "Baby check",
                                "Baby check & maternal check",
                                "Baby check & infant vaccination",
                                "All 3 services")) 

# calculate percentages of each combo in each ethnic group
eth_denom_main_allthree <- babycheck_mainanalysis_allthree %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  ungroup() %>% 
  rename(N = n)
eth_main_allthree <- babycheck_mainanalysis_allthree %>% 
  group_by(M_eth18_2011, allthree) %>% 
  count() %>% 
  ungroup()
eth_main_allthree <- eth_main_allthree %>% 
  left_join(eth_denom_main_allthree, by = c("M_eth18_2011")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(percent_unrounded = n/N*100) %>% 
  mutate(lowerci = prop-1.96*sqrt((prop*(1-prop))/N)) %>%
  mutate(upperci = prop+1.96*sqrt((prop*(1-prop))/N)) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1))

write_csv(eth_main_allthree, file = "~/eth_main_allthree.csv")

# stacked bar chart by ethnic group
eth_main_allthree_stacked <- ggplot(data = eth_main_allthree, aes(fill=allthree, x = percent_unrounded, y = M_eth18_2011)) +
  geom_bar(stat = "identity") + 
  guides(fill=guide_legend(title = "Services received on the same day")) +
  theme(text = element_text(size =10)) +
  theme_gray(base_size = 20) + # scales size of the plot relative to the text
  labs(x = "Percent", y = "Maternal ethnicity", title = " ") +
  scale_fill_manual(values=c("lightskyblue1","skyblue2", "dodgerblue3", "dodgerblue4")) +
  scale_y_discrete(limits = rev)  # flip y axis
eth_main_allthree_stacked

ggsave(width = 18, height = 8, dpi = 450, "~/eth_main_allthree.jpg")

# Baby check + maternal check + infant vax: sens --------------------------------

# Repeat for babycheck_sensanalysis as necessary 
