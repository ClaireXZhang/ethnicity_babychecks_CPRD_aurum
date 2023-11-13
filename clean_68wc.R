#--------------------------------------------------------------------------
# SET UP                          
#--------------------------------------------------------------------------
# Description -------------------------------------------------------------

# Cleaning of 6-8 week baby checks and maternal checks in CPRD Aurum

# Set working directory ---------------------------------------------------

setwd('K:/')

# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse, lubridate)

# Disable scientific notation ---------------------------------------------

options(scipen = 999) 

#--------------------------------------------------------------------------
# IDENTIFY 6-8 WEEK BABY CHECKS                          
#--------------------------------------------------------------------------
# Load --------------------------------------------------------------------

# load child observations file
load(file="~/obs_C_clean.Rdata")
load(file="~/patids.Rdata")

# filter by those in the patid_excl list (full cohort after all exclusions) 
obs_C_clean <- obs_C_clean %>% 
  filter(C_patid %in% patids$C_patid) 

# 6-8 week baby check codelist 
babycheck_codelist <- read_tsv('~/six_week_baby_check_medcodes_final.txt', 
                                  col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) 

# 6-8 week baby check refusal code list
babycheck_refusal_codelist <- read_tsv('~/six_week_baby_check_refusal.txt', 
                              col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) 

# Identify baby checks in obs file ----------------------------------------

# identify obs with 6 week baby check codes
obs_C_babycheck <- obs_C_clean %>% 
  inner_join(babycheck_codelist, by = "medcodeid") %>% 
  select(C_patid, obsdate, category, checklist, consid) %>% 
  mutate(babycheck = 1)

# check proportion of obs that are checklist codes
total_obs <- nrow(obs_C_babycheck)
n_checklist_obs <- sum(!is.na(obs_C_babycheck$checklist))
percent_checklist_obs <- round(n_checklist_obs/total_obs*100, 2) 

# remove checklist variable, not needed from hereon
obs_C_babycheck <- obs_C_babycheck %>% 
  select(-checklist)

# save
save(obs_C_babycheck, file="~/obs_C_babycheck.Rdata")

# Counts & exclusions ------------------------------------------------

# count obs with no obsdate
nodate_obs <- sum(is.na(obs_C_babycheck$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# identify and count refusal obs
babycheck_refusal_codelist <- babycheck_refusal_codelist %>% 
  select(medcodeid) 
obs_sixweekrefusal <- obs_C_clean %>% 
  inner_join(babycheck_refusal_codelist, by = "medcodeid") %>% 
  select(C_patid, obsdate)
obs_sixweekrefusal <- obs_sixweekrefusal %>% 
  mutate(combine = paste0(C_patid,":",obsdate))
obs_C_babycheck <- obs_C_babycheck %>% 
  mutate(combine = paste0(C_patid,":",obsdate))
n_refusal_obs <- sum(obs_C_babycheck$combine %in% obs_sixweekrefusal$combine) 
percent_refusal_obs <- round(n_refusal_obs/total_obs*100, 2)

# drop obs with no obsdate and refusal obs 
obs_C_babycheck <- obs_C_babycheck %>% 
  filter(!is.na(obs_C_babycheck$obsdate)) %>% 
  filter(!(combine %in% obs_sixweekrefusal$combine)) %>% 
  select(-c(combine))

# count obs after exclusions
total_obs_excl <- nrow(obs_C_babycheck)

# Joint to patid list, count obs occurring outside of 4-12 week window -----------------------------------------------------

# load C_patid list - TO UPDATE
load(file = "~/patids.Rdata")

babypatids <- patids %>% 
  select(C_patid, deldate) 

# join to C_patid list
C_babycheck <- babypatids %>% 
  left_join(obs_C_babycheck, by = "C_patid") %>% 
  mutate(babycheck = replace_na(babycheck, 0))

# Check that those without babycheck codes are distinct patients, and put them into a separate df to join later
C_no_babycheck <- C_babycheck %>% 
  filter(babycheck == 0)
n_distinct(C_no_babycheck$C_patid) == nrow(C_no_babycheck)

# Put the rest (patients with babycheck codes) into separate df
C_yes_babycheck <- C_babycheck %>% 
  filter(babycheck == 1)

# check that these split dfs still have the same number of unique patients altogether
n_distinct(C_yes_babycheck$C_patid) + n_distinct(C_no_babycheck$C_patid) == n_distinct(C_babycheck$C_patid)

# count obs with dates outside of 4-12 weeks after birth (for consistency, denominator used is the one before exclusion)
# drop these from C_yes_babycheck and move these into a separate df, select the first row to get one row per patient and join with the patients with no babycheck codes
n_dates_outside_4to12wk <- sum(C_yes_babycheck$obsdate < (C_yes_babycheck$deldate %m+% weeks(4)) | C_yes_babycheck$obsdate >= (C_yes_babycheck$deldate%m+% weeks(13)))
percent_dates_outside_4to12wk <- round(n_dates_outside_4to12wk/total_obs*100, 2) 

C_babycheck_dates_outside_4to12wk <- C_yes_babycheck %>%
  filter(obsdate < (deldate %m+% weeks(4)) | obsdate >= (deldate %m+% weeks(13))) %>%
  mutate(babycheck = 0) %>%
  group_by(C_patid) %>%
  slice(1)

C_no_babycheck <- rbind(C_no_babycheck, C_babycheck_dates_outside_4to12wk) 
n_distinct(C_no_babycheck$C_patid) == nrow(C_no_babycheck) # each row is a unique patient

C_yes_babycheck_dates <- C_yes_babycheck %>%
  filter(obsdate >= (deldate %m+% weeks(4))) %>%
  filter(obsdate < (deldate %m+% weeks(13)))

# Find obs with 1 unique code and split these into separate df
C_yes_babycheck_unique <- C_yes_babycheck_dates %>% 
  group_by(C_patid) %>%
  count() %>% 
  filter(n == 1) %>% 
  select(-n) %>% 
  left_join(C_yes_babycheck_dates, by = "C_patid") # only have 1 code
n_distinct(C_yes_babycheck_unique$C_patid) == nrow(C_yes_babycheck_unique) # each row is a unique patient

# Prioritise by level of certainty to cut down to 1 record per person --------

# identify patients with more than 1 code
# create priority variable and only keep the row with the highest priority for each patient 
C_yes_babycheck_unique2 <- C_yes_babycheck_dates %>%
  filter(!(C_patid %in% C_yes_babycheck_unique$C_patid)) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, category == "dev_check", 2)) %>% 
  mutate(priority = replace(priority, category == "6wk", 3)) %>% 
  mutate(priority = replace(priority, category == "health_visitor", 4)) %>% 
  group_by(C_patid) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority)
n_distinct(C_yes_babycheck_unique2$C_patid) == nrow(C_yes_babycheck_unique2) # each row is a unique patient

# remove patients from C_no_babycheck if they appear in the other 2 dfs with unique babycheck codes
# set other fields to NA
C_no_babycheck <- C_no_babycheck %>% 
  filter(!(C_patid %in% C_yes_babycheck_unique$C_patid)) %>% 
  filter(!(C_patid %in% C_yes_babycheck_unique2$C_patid)) %>% 
  mutate(obsdate = NA) %>%  
  mutate(consid = NA) %>%
  mutate(category = NA)
  
# join this to the previous dfs of patients with unique codes to begin with and those with no codes
C_babycheck_clean <- rbind(C_yes_babycheck_unique2, C_yes_babycheck_unique, C_no_babycheck)
n_distinct(C_babycheck_clean$C_patid) == nrow(C_babycheck_clean) # each row is a unique patient

# check this final df has the same number of patients as the patid list
nrow(patids) == nrow(C_babycheck_clean) 

# save
save(C_babycheck_clean, file="~/C_babycheck_clean.Rdata")

# Find 6WCs that are linked to consultations within 4-12 weeks ------------

# load
load(file="~/cons_direct_C_clean.Rdata")

# to get one row per consultation, group by consdate and consid, prioritising staff in order of GP, nurse, associate, community/hv
cons_direct_C_clean_dedup <- cons_direct_C_clean %>%
  mutate(job_priority = 1) %>% 
  mutate(job_priority = replace(job_priority, job == "Nurse", 2)) %>% 
  mutate(job_priority = replace(job_priority, job == "Associate/practitioner", 3)) %>% 
  mutate(job_priority = replace(job_priority, job == "HV/midwife/community", 4))

cons_direct_C_clean_dedup <- cons_direct_C_clean_dedup %>%
  group_by(C_patid, consdate, consid) %>%
  slice(which.min(job_priority)) %>%
  mutate(cons = 1)

n_distinct(cons_direct_C_clean_dedup$consid) == nrow(cons_direct_C_clean_dedup) # each row is a distinct consid

# join to 68wc data
C_babycheck_yes_clean <- C_babycheck_clean %>% 
  filter(babycheck == 1) %>% 
  left_join(cons_direct_C_clean_dedup, by = c("C_patid", "consid", "obsdate" = "consdate"))

n_distinct(C_babycheck_yes_clean$C_patid) == nrow(C_babycheck_yes_clean)

# deselect unneeded variables
C_babycheck_yes_clean <- C_babycheck_yes_clean %>% 
  select(-c(job_priority, job_detail, cons_description, staffid, C_pracid))

# check how many consultations were linked
sum(C_babycheck_yes_clean$cons == 1, na.rm=T) 

# create NA values for cons-related variables for people who had no evidence of direct consultations
C_babycheck_no_clean <- C_babycheck_clean %>% 
  filter(babycheck == 0) 
C_babycheck_no_clean$job <- as.character(NA)
C_babycheck_no_clean$f2f <- as.character(NA)
C_babycheck_no_clean$cons <- as.numeric(NA)

# join people who had and didn't have direct consultations
C_babycheck_cons_clean <- rbind(C_babycheck_yes_clean, C_babycheck_no_clean) %>% 
  mutate(cons = replace_na(cons, 0)) 

# check this new df has the same number of patients as the original one
n_distinct(C_babycheck_cons_clean$C_patid) == n_distinct(C_babycheck_clean$C_patid) 

# save
save(C_babycheck_cons_clean, file="~/C_babycheck_cons_clean.Rdata")

# Save counts and checks  ---------------------------------------------------

total_obs <- as.data.frame(total_obs)
n_checklist_obs <- as.data.frame(n_checklist_obs)
percent_checklist_obs <- as.data.frame(percent_checklist_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
n_refusal_obs <- as.data.frame(n_refusal_obs)
percent_refusal_obs <- as.data.frame(percent_refusal_obs)
n_dates_outside_4to12wk <- as.data.frame(n_dates_outside_4to12wk)
percent_dates_outside_4to12wk <- as.data.frame(percent_dates_outside_4to12wk)
total_obs_excl <- as.data.frame(total_obs_excl)

babycheck_cleaning_counts <- cbind(total_obs, n_checklist_obs,percent_checklist_obs,
                                   nodate_obs, percent_nodate_obs,
                                   n_refusal_obs, percent_refusal_obs,
                                   n_dates_outside_4to12wk, percent_dates_outside_4to12wk,
                                   total_obs_excl)

write_csv(babycheck_cleaning_counts, file = "~/babycheck_cleaning_counts.csv")

rm(list=ls())

#--------------------------------------------------------------------------
# IDENTIFY 6-8 WEEK MATERNAL CHECKS                          
#--------------------------------------------------------------------------
# Load --------------------------------------------------------------------

# load child observations and drug issue file
load(file="~/obs_M_clean.Rdata")
load(file="~/patids.Rdata")

# filter by those in the patid_excl list  
obs_M_clean <- obs_M_clean %>% 
  filter(patid %in% patids$M_patid) 

# 6-8 week maternal check codelist 
maternalcheck_codelist <- read_tsv('~/maternal_check_codes_medcodes_all.txt', 
                                              col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) 

# exclude irrelevant baby check codes (kept in the txt file for now just in case we want to use them later)
maternalcheck_codelist <- maternalcheck_codelist %>%
  select(medcodeid, category)

# Identify maternal checks in obs file ----------------------------------------

# identify obs with 6 week baby check codes
obs_M_maternalcheck <- obs_M_clean %>% 
  inner_join(maternalcheck_codelist, by = "medcodeid") %>% 
  select(M_patid, obsdate, category, consid) %>% 
  mutate(maternalcheck = 1)

# save
save(obs_M_maternalcheck, file="~/obs_M_maternalcheck.Rdata")

# Counts & exclusions ------------------------------------------------

total_obs <- nrow(obs_M_maternalcheck)

# count obs with no obsdate
nodate_obs <- sum(is.na(obs_M_maternalcheck$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# drop obs with no obsdate 
obs_M_maternalcheck <- obs_M_maternalcheck %>% 
  filter(!is.na(obsdate)) 

# count obs after exclusions
total_obs_excl <- nrow(obs_M_maternalcheck)

# Joint to patid list, count obs occurring outside of 4-12 week window -----------------------------------------------------

# load M_patid list 
load(file = "~/patids.Rdata")

# join to M_patid list
M_maternalcheck <- patids %>% 
  left_join(obs_M_maternalcheck, by = "M_patid") %>% 
  mutate(maternalcheck = replace_na(maternalcheck, 0))

# Check that those without maternalcheck codes are distinct patients, and put them into a separate df to join later
M_no_maternalcheck <- M_maternalcheck %>% 
  filter(maternalcheck == 0)
n_distinct(M_no_maternalcheck$C_patid) == nrow(M_no_maternalcheck) # one unique child per row (not unique mothers because some mums have multiple children)

# Put the rest (patients with maternalcheck codes) into separate df
M_yes_maternalcheck <- M_maternalcheck %>% 
  filter(maternalcheck == 1)

# check that these split dfs still have the same number of unique patients altogether
n_distinct(M_yes_maternalcheck$M_patid) + n_distinct(M_no_maternalcheck$M_patid) == n_distinct(M_maternalcheck$M_patid)
n_distinct(M_yes_maternalcheck$C_patid) + n_distinct(M_no_maternalcheck$C_patid) == n_distinct(M_maternalcheck$C_patid)

# count obs with dates outside of 4-12 weeks after birth (for consistency, denominator used is the one before exclusion)
# drop these from M_yes_maternalcheck and move these into a separate df, select the first row to get one row per patient and join with the patients with no maternalcheck codes
n_dates_outside_4to12wk <- sum(M_yes_maternalcheck$obsdate < (M_yes_maternalcheck$deldate %m+% weeks(4)) | M_yes_maternalcheck$obsdate >= (M_yes_maternalcheck$deldate%m+% weeks(13)))
percent_dates_outside_4to12wk <- round(n_dates_outside_4to12wk/total_obs*100, 2) # we haven't differentiated between 6-8 week maternal check for a child and their sibling at this point

M_maternalcheck_dates_outside_4to12wk <- M_yes_maternalcheck %>%
  filter(obsdate < (deldate %m+% weeks(4)) | obsdate >= (deldate %m+% weeks(13))) %>%
  mutate(maternalcheck = 0) %>%
  group_by(C_patid) %>%
  slice(1)

M_no_maternalcheck <- rbind(M_no_maternalcheck, M_maternalcheck_dates_outside_4to12wk) 
n_distinct(M_no_maternalcheck$C_patid) == nrow(M_no_maternalcheck) # each row is a unique baby

M_yes_maternalcheck_dates <- M_yes_maternalcheck %>%
  filter(obsdate >= (deldate %m+% weeks(4))) %>%
  filter(obsdate < (deldate %m+% weeks(13)))

# Find obs with 1 unique code and split these into separate df
M_yes_maternalcheck_unique <- M_yes_maternalcheck_dates %>% 
  group_by(C_patid) %>%
  count() %>% 
  filter(n == 1) %>% 
  select(-n) %>% 
  left_join(M_yes_maternalcheck_dates, by = "C_patid") # only have 1 code
n_distinct(M_yes_maternalcheck_unique$C_patid) == nrow(M_yes_maternalcheck_unique) # each row is a unique baby

# Prioritise by level of certainty to cut down to 1 record per person --------

# identify patients with more than 1 code
# create priority variable and only keep the row with the highest priority for each patient
M_yes_maternalcheck_unique2 <- M_yes_maternalcheck_dates %>%
  filter(!(C_patid %in% M_yes_maternalcheck_unique$C_patid)) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, category == "postnatal", 2)) %>% 
  mutate(priority = replace(priority, category == "6wk", 3)) %>% 
  mutate(priority = replace(priority, category == "health_visitor", 4)) %>% 
  group_by(C_patid) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority)
n_distinct(M_yes_maternalcheck_unique2$C_patid) == nrow(M_yes_maternalcheck_unique2) # each row is a unique baby

# remove patients from M_no_maternalcheck if they appear in the other 2 dfs with unique maternalcheck codes
# set other fields to NA
M_no_maternalcheck <- M_no_maternalcheck %>% 
  filter(!(C_patid %in% M_yes_maternalcheck_unique$C_patid)) %>% 
  filter(!(C_patid %in% M_yes_maternalcheck_unique2$C_patid)) %>% 
  mutate(obsdate = NA) %>%  
  mutate(consid = NA) %>%
  mutate(category = NA)

# join this to the previous dfs of patients with unique codes to begin with and those with no codes
M_maternalcheck_clean <- rbind(M_yes_maternalcheck_unique2, M_yes_maternalcheck_unique, M_no_maternalcheck)
n_distinct(M_maternalcheck_clean$C_patid) == nrow(M_maternalcheck_clean) # each row is a unique patient

# check this final df has the same number of patients as the patid list
nrow(patids) == nrow(M_maternalcheck_clean) 

# save
save(M_maternalcheck_clean, file="~/M_maternalcheck_clean.Rdata")

# Find 6WCs that are linked to consultations within 4-12 weeks ------------

# load
load(file="~/cons_direct_M_clean.Rdata")

# to get one row per consultation, group by consdate and consid, prioritising staff in order of GP, nurse, associate, community/hv
cons_direct_M_clean_dedup <- cons_direct_M_clean %>%
  mutate(job_priority = 1) %>% 
  mutate(job_priority = replace(job_priority, job == "Nurse", 2)) %>% 
  mutate(job_priority = replace(job_priority, job == "Associate/practitioner", 3)) %>% 
  mutate(job_priority = replace(job_priority, job == "HV/midwife/community", 4))

cons_direct_M_clean_dedup <- cons_direct_M_clean_dedup %>%
  group_by(M_patid, consdate, consid) %>%
  slice(which.min(job_priority)) %>%
  mutate(cons = 1)

n_distinct(cons_direct_M_clean_dedup$consid) == nrow(cons_direct_M_clean_dedup) # each row is a distinct consid

# join to 68wc data
M_maternalcheck_yes_clean <- M_maternalcheck_clean %>% 
  filter(maternalcheck == 1) %>% 
  left_join(cons_direct_M_clean_dedup, by = c("M_patid", "consid", "obsdate" = "consdate"))

n_distinct(M_maternalcheck_yes_clean$C_patid) == nrow(M_maternalcheck_yes_clean)

# deselect unneeded variables
M_maternalcheck_yes_clean <- M_maternalcheck_yes_clean %>% 
  select(-c(job_priority, job_detail, cons_description, staffid, M_pracid))

# check how many consultations were linked
sum(M_maternalcheck_yes_clean$cons == 1, na.rm=T) 

# create NA values for cons-related variables for people who had no evidence of direct consultations
M_maternalcheck_no_clean <- M_maternalcheck_clean %>% 
  filter(maternalcheck == 0) 
M_maternalcheck_no_clean$job <- as.character(NA)
M_maternalcheck_no_clean$f2f <- as.character(NA)
M_maternalcheck_no_clean$cons <- as.numeric(NA)

# join people who had and didn't have direct consultations
M_maternalcheck_cons_clean <- rbind(M_maternalcheck_yes_clean, M_maternalcheck_no_clean) %>% 
  mutate(cons = replace_na(cons, 0)) 

# check this new df has the same number of patients as the original one
n_distinct(M_maternalcheck_cons_clean$M_patid) == n_distinct(M_maternalcheck_clean$M_patid) 
n_distinct(M_maternalcheck_cons_clean$C_patid) == n_distinct(M_maternalcheck_clean$C_patid) 

# save
save(M_maternalcheck_cons_clean, file="~/M_maternalcheck_cons_clean.Rdata")

# Apply limits for mum's GP registration dates --------

# load
load(file="~/aurum_CM_clean.Rdata")
load(file="~/C_babycheck_cons_clean.Rdata")

# select needed variables
fu_checks <- aurum_CM_clean %>% 
  select(c(M_patid, C_patid, M_regstartdate, M_pc_fu_end))

# count how many mothers were not registered from 4 weeks to the end of 12 weeks
M_maternalcheck_cons_clean <- M_maternalcheck_cons_clean %>% 
  left_join(fu_checks, by = c("M_patid", "C_patid")) 
  
total_child <- n_distinct(C_babycheck_cons_clean$C_patid)
not_registered_4to12wk <- sum(!(M_maternalcheck_cons_clean$M_regstartdate <= M_maternalcheck_cons_clean$deldate %m+% weeks(4)& 
                                      M_maternalcheck_cons_clean$M_pc_fu_end >= M_maternalcheck_cons_clean$deldate %m+% weeks(13)),
                                  na.rm = T)
percent_not_registered_4to12wk <- round(not_registered_4to12wk/total_child*100, 2) # denominator is number of children affected by this

# identify mothers who were not registered from 4 weeks to the end of 12 weeks
# set these records to maternalcheck_4to12reg = 0 
M_maternalcheck_cons_clean <- M_maternalcheck_cons_clean %>% 
  mutate(maternalcheck_4to12reg = ifelse(M_pc_fu_end >= deldate %m+% weeks(13) & M_regstartdate <= deldate %m+% weeks(4), maternalcheck, 0)) 

# count how many mothers were not registered by 12 months or had follow-up end before end of 12 weeks
not_registered_by12m_fuendbefore12wk <- sum(!(M_maternalcheck_cons_clean$M_regstartdate <= M_maternalcheck_cons_clean$deldate %m+% years(1)& 
                                      M_maternalcheck_cons_clean$M_pc_fu_end >= M_maternalcheck_cons_clean$deldate %m+% weeks(13)),
                                  na.rm = T)
percent_not_registered_by12m_fuendbefore12wk <- round(not_registered_by12m_fuendbefore12wk/total_child*100, 2) # denominator is number of children affected by this

# identify mums whose follow-up didn't end before end of 12 weeks and reg start can be anytime up to 12 months
# set these records to maternalcheck = 0 
M_maternalcheck_cons_clean <- M_maternalcheck_cons_clean %>%
  mutate(maternalcheck_12m = ifelse(M_pc_fu_end >= deldate %m+% weeks(13) & M_regstartdate <= deldate %m+% years(1), maternalcheck, 0)) 

# save
save(M_maternalcheck_cons_clean, file="~/M_maternalcheck_cons_clean.Rdata")

# Save counts and checks ---------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
n_dates_outside_4to12wk <- as.data.frame(n_dates_outside_4to12wk)
percent_dates_outside_4to12wk <- as.data.frame(percent_dates_outside_4to12wk)
total_obs_excl <- as.data.frame(total_obs_excl)
not_registered_birthto12wk <- as.data.frame(not_registered_birthto12wk)
percent_not_registered_birthto12wk <- as.data.frame(percent_not_registered_birthto12wk)
not_registered_by12m_fuendbefore12wk <- as.data.frame(not_registered_by12m_fuendbefore12wk)
percent_not_registered_by12m_fuendbefore12wk <- as.data.frame(percent_not_registered_by12m_fuendbefore12wk)

maternalcheck_cleaning_counts <- cbind(total_obs,
                                   nodate_obs, percent_nodate_obs,
                                   n_dates_outside_4to12wk, percent_dates_outside_4to12wk,
                                   total_obs_excl,
                                   not_registered_birthto12wk, percent_not_registered_birthto12wk,
                                   not_registered_by12m_fuendbefore12wk, percent_not_registered_by12m_fuendbefore12wk)

write_csv(maternalcheck_cleaning_counts, file = "~/maternalcheck_cleaning_counts.csv")

rm(list=ls())


#--------------------------------------------------------------------------
# SERVICE DELIVERY                           
#--------------------------------------------------------------------------
# Derive baby check - maternal check variable  -----------------------------

# join baby check and maternal check files
load(file="~/C_babycheck_cons_clean.Rdata")
load(file="~/M_maternalcheck_cons_clean.Rdata")

M_maternalcheck_join <- M_maternalcheck_cons_clean %>% 
  select(M_patid, maternalcheck, maternalcheck_4to12reg, maternalcheck_12m, C_patid, obsdate) %>% 
  rename(M_obsdate = obsdate)

# derive babycheck_maternalcheck for mums regardless of registration (1 = maternal check before baby check, 2 = both on the same date, 3 = maternal check after baby check, 4 = no maternal check)
C_babycheck_outcomes <- C_babycheck_cons_clean %>% 
  left_join(M_maternalcheck_join, by = c("C_patid")) %>% 
  mutate(babycheck_maternalcheck = 4) %>% 
  mutate(babycheck_maternalcheck = replace(babycheck_maternalcheck, M_obsdate < obsdate, 1)) %>% 
  mutate(babycheck_maternalcheck = replace(babycheck_maternalcheck, M_obsdate == obsdate, 2)) %>% 
  mutate(babycheck_maternalcheck = replace(babycheck_maternalcheck, M_obsdate > obsdate, 3)) %>% 
  mutate(babycheck_maternalcheck = replace(babycheck_maternalcheck, babycheck == 0, NA))

# derive babycheck_maternalcheck for mums registered between 4-12 weeks 
C_babycheck_outcomes <- C_babycheck_outcomes %>% 
  mutate(babycheck_maternalcheck_4to12reg = ifelse(maternalcheck_4to12reg == 0 & babycheck == 1, 4, babycheck_maternalcheck))
  
# derive babycheck_maternalcheck for mums whose follow-up didn't end before 12 weeks, but reg start can be anytime up to 12 months
C_babycheck_outcomes <- C_babycheck_outcomes %>%
  mutate(babycheck_maternalcheck_12m = ifelse(maternalcheck_12m == 0 & babycheck == 1, 4, babycheck_maternalcheck))

# save
save(C_babycheck_outcomes, file="~/C_babycheck_outcomes.Rdata")

# Derive baby check - vaccination variable  -------------

# load baby check analysis file and all vaccination files for vaccinations expected in early infancy
# see vaccinations GitHub repo (link in README file)
load(file="~/C_babycheck_outcomes.Rdata")
load(file="~/obs_rotavirus_full.Rdata")
load(file="~/obs_pneumococcal_full.Rdata")
load(file="~/obs_menb_full.Rdata")
load(file="~/obs_menc_full.Rdata")
load(file="~/obs_6in1_full.Rdata")
load(file="~/obs_5in1_full.Rdata")

# select needed variables
obs_rotavirus_full <- obs_rotavirus_full %>% dplyr::select(C_patid, obsdate)
obs_pneumococcal_full <- obs_pneumococcal_full %>% dplyr::select(C_patid, obsdate)
obs_menb_full <- obs_menb_full %>% dplyr::select(C_patid, obsdate)
obs_menc_full <- obs_menc_full %>% dplyr::select(C_patid, obsdate)
obs_6in1_full <- obs_6in1_full %>% dplyr::select(C_patid, obsdate)
obs_5in1_full <- obs_5in1_full %>% dplyr::select(C_patid, obsdate)

# join them 
obs_infant_vax <- rbind(obs_rotavirus_full, obs_pneumococcal_full,
                        obs_menb_full, obs_menc_full,
                        obs_5in1_full, obs_6in1_full)
  
# drop those with missing obsdate, and deduplicate to one vax per obsdate
obs_infant_vax <- obs_infant_vax %>% 
  filter(!is.na(obsdate)) %>% 
  distinct(C_patid, obsdate, .keep_all = TRUE) %>% 
  mutate(babycheck_vax = 1)

# derive babycheck_vax (0 = 6-8 week baby check in isolation, 1 = together with any vaccination with on the same date)
C_babycheck_outcomes <- C_babycheck_outcomes %>% 
  left_join(obs_infant_vax, by = c("C_patid", "obsdate")) %>% 
  mutate(babycheck_vax = replace(babycheck_vax, babycheck == 1 & is.na(babycheck_vax), 0))

# derive variable indicating that infant vax occurred at any point during 4-12 week period
obs_infant_vax <- obs_infant_vax %>% 
  mutate(infantvax = 1) %>% 
  rename(obsdatevax = obsdate) %>% 
  select(C_patid, infantvax, obsdatevax)
C_babycheck_outcomes <- C_babycheck_outcomes %>% 
  left_join(obs_infant_vax, by = "C_patid") %>% 
  mutate(infantvax = ifelse(obsdatevax < deldate %m+% weeks(4) | obsdatevax >= deldate %m+% weeks(13), 0, infantvax)) %>% 
  select(-obsdatevax) %>% 
  mutate(infantvax = replace_na(infantvax, 0)) %>% 
  group_by(C_patid) %>% 
  slice(which.max(infantvax)) %>% 
  ungroup()

# save
save(C_babycheck_outcomes, file="~/C_babycheck_outcomes.Rdata")

table(C_babycheck_outcomes$babycheck_vax)
table(C_babycheck_outcomes$infantvax)

#--------------------------------------------------------------------------
# JOIN ALL REQUIRED VARIABLES INTO FULL ANALYSIS FILE                
#--------------------------------------------------------------------------

# Load all files ----------------------------------------------------------

# see ethnicity & vaccinations GitHub repo (link in README file)
load(file="~/CM_ethnicity_clean.Rdata")
load(file="~/aurum_CM_clean.Rdata")
load(file="~/delvar_all_clean.Rdata")
load(file="~/arealevel_CM_clean.Rdata")
load(file="~/CM_hospitalised_babycheck.Rdata")
load(file="~/C_babycheck_outcomes.Rdata")

# select needed variables (no duplication of variables across files)
aurum <- aurum_CM_clean %>% 
  select(C_patid, C_pracid, C_gender, C_region, dob_deldate, C_regstartdate, C_pc_fu_end, M_patid, matage_aurum_cat, M_regstartdate)
ethnicity <- CM_ethnicity_clean %>% 
  select(C_patid, C_eth18_2011, C_eth6_2011, C_eth5_2011, M_eth18_2011, M_eth6_2011, M_eth5_2011, CM_eth_discord)
delvar <- delvar_all_clean %>% 
  select(C_patid, gest_cat, bw_cat, first_time_mum, delmeth_cat)
arealevel <- arealevel_CM_clean %>% 
  select(C_patid, imd, ruc)
babycheck <- C_babycheck_outcomes %>% 
  select(C_patid, babycheck, category, job, f2f, obsdate,
         maternalcheck, maternalcheck_4to12reg, maternalcheck_12m,
         babycheck_maternalcheck, babycheck_maternalcheck_4to12reg, babycheck_maternalcheck_12m, babycheck_vax)

# Join all variables and create time period variable ---------------------------------------------

# join all variables
babycheck_analysis <- babycheck %>% 
  left_join(aurum, by = "C_patid") %>% 
  left_join(ethnicity, by = "C_patid") %>% 
  left_join(arealevel, by = "C_patid") %>% 
  left_join(delvar, by = "C_patid") %>% 
  left_join(CM_hospitalised_babycheck, by = "C_patid") 
  
# check 1 row per child
n_distinct(babycheck_analysis$C_patid) == nrow(babycheck_analysis)

# check number of children in analysis file = outcome file 
n_distinct(babycheck_analysis$C_patid) == n_distinct(C_babycheck_outcomes$C_patid)

# create time period variable (financial year in which baby was born)
babycheck_analysis <- babycheck_analysis %>% 
  mutate(financialyear = NA) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate < "2006-04-01", "2005-06")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2006-04-01" & dob_deldate < "2007-04-01", "2006-07")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2007-04-01" & dob_deldate < "2008-04-01", "2007-08")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2008-04-01" & dob_deldate < "2009-04-01", "2008-09")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2009-04-01" & dob_deldate < "2010-04-01", "2009-10")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2010-04-01" & dob_deldate < "2011-04-01", "2010-11")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2011-04-01" & dob_deldate < "2012-04-01", "2011-12")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2012-04-01" & dob_deldate < "2013-04-01", "2012-13")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2013-04-01" & dob_deldate < "2014-04-01", "2013-14")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2014-04-01" & dob_deldate < "2015-04-01", "2014-15")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2015-04-01" & dob_deldate < "2016-04-01", "2015-16")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2016-04-01" & dob_deldate < "2017-04-01", "2016-17")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2017-04-01" & dob_deldate < "2018-04-01", "2017-18")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2018-04-01" & dob_deldate < "2019-04-01", "2018-19")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2019-04-01" & dob_deldate < "2020-04-01", "2019-20")) %>% 
  mutate(financialyear = replace(financialyear, dob_deldate >= "2020-04-01", "2020-21"))  
  
# rename variables for ease of labelling during analysis and set to correct variable type
babycheck_analysis <- babycheck_analysis %>% 
  mutate(financialyear = as.factor(financialyear)) %>% 
  mutate(category = as.factor(category)) %>% 
  rename(certainty = category) %>% 
  mutate(f2f = as.numeric(f2f)) %>% 
  rename(hp = job) %>% 
  rename(sex = C_gender) %>% 
  rename(region = C_region) %>% 
  rename(dob = dob_deldate) %>% 
  mutate(matage_aurum_cat = as.factor(matage_aurum_cat)) %>% 
  rename(matage = matage_aurum_cat) %>% 
  mutate(gest_cat = as.factor(gest_cat)) %>% 
  rename(gest = gest_cat) %>% 
  mutate(bw_cat = as.factor(bw_cat)) %>% 
  rename(bw = bw_cat) %>% 
  mutate(delmeth_cat = as.factor(delmeth_cat)) %>% 
  rename(delmeth = delmeth_cat) 
  
summary(is.na(babycheck_analysis))

# relevel factors
babycheck_analysis <- babycheck_analysis %>% 
  mutate(financialyear = fct_relevel(financialyear, "2005-06", "2006-07","2007-08","2008-09","2009-10","2010-11",
                                     "2011-12","2012-13","2013-14","2014-15","2015-16","2016-17","2017-18","2018-19",
                                     "2019-20","2020-21")) %>%
  mutate(certainty = fct_relevel(certainty,"6wk_baby_check", "6wk", "dev_check", "health_visitor")) %>%
  mutate(hp = fct_relevel(hp, "GP/doctor", "Nurse", "Associate/practitioner", "HV/midwife/community")) %>%
  mutate(sex = fct_relevel(sex, "M", "F", "I")) %>%
  mutate(region = fct_relevel(region, "London", "East Midlands" ,"East of England","North East" ,
                              "North West","South East" ,"South West" ,"West Midlands","Yorkshire and The Humber")) %>%
  mutate(matage = fct_relevel(matage, "Under 20","20-24", "25-29", "30-34", "35-39", "40 and older")) %>% 
  mutate(gest = fct_relevel(gest,"extremely preterm (<28 weeks)", "very preterm (28-31 weeks)",  
                            "moderate to late preterm (32-36 weeks)", "term (37-42 weeks)")) %>%
  mutate(bw = fct_relevel(bw,"very low", "low", "normal", "high")) %>%
  mutate(delmeth = fct_relevel(delmeth,"spontaneous vaginal", "instrumental or other", 
                               "elective caesarean", "emergency caesarean")) 
  
# change numerical variables to factors (except babycheck outcome - not allowed for poisson function)
babycheck_analysis <- babycheck_analysis %>% 
  mutate(across(c(where(is.numeric), -babycheck), as.factor))
  
# save
save(babycheck_analysis, file="~/babycheck_analysis.Rdata")

# Create aggregate birth variables for modelling --------------------------

# preterm vs term
babycheck_analysis <- babycheck_analysis %>% 
  mutate(preterm = 1) %>% 
  mutate(preterm = replace(preterm, gest == "term (37-42 weeks)", 0))

# spontaneous vaginal vs other modes of birth
babycheck_analysis <- babycheck_analysis %>% 
  mutate(modebirth = 1) %>% 
  mutate(modebirth = replace(modebirth, delmeth == "spontaneous vaginal", 0))

# save
save(babycheck_analysis, file="~/babycheck_analysis.Rdata")

#--------------------------------------------------------------------------
# FLAGS & EXCLUSIONS FOR 6-8 WEEK CHECKS OUTCOMES                       
#--------------------------------------------------------------------------
# Follow-up and registration start dates -----------------------

# load
load(file="~/babycheck_analysis.Rdata")

# how many baby checks were conducted 4-5 weeks (to assess age cutoff for registration)
total_68wc <- sum(babycheck_analysis$babycheck == 1)
n_obsdate_45wk <- sum(babycheck_analysis$obsdate >= babycheck_analysis$dob %m+% weeks(4) & 
                        babycheck_analysis$obsdate < babycheck_analysis$dob %m+% weeks(6), na.rm = T)
percent_obsdate_45wk <- round(n_obsdate_45wk/total_68wc*100, 2) 

# count how many children have end of follow-up before end of 12 weeks
total_child <- nrow(babycheck_analysis)
fuend_before12wk <- sum(babycheck_analysis$C_pc_fu_end < (babycheck_analysis$dob %m+% weeks(13)))
percent_fuend_before12wk <- round(fuend_before12wk/total_child*100, 2) 

# Count & make exclusions ---------------------------------------------------------

# indeterminate sex
indeterminate_sex <- sum(babycheck_analysis$sex == "I")
percent_indeterminate_sex <- round(indeterminate_sex/total_child*100, 3)
sourcepop <- babycheck_analysis %>% 
  filter(sex != "I") 

# multiple birth
load(file="~/delvar_all_clean.Rdata")
multiple_birth <- sum(delvar_all_clean$multiple == 1)
percent_multiple_birth <- round(multiple_birth/total_child*100, 3)
singleton <- delvar_all_clean %>% 
  filter(multiple == 0) %>% 
  select(C_patid)
sourcepop <- sourcepop %>% 
  filter(C_patid %in% singleton$C_patid)

# missing region
missing_region <- sum(is.na(babycheck_analysis$region))
percent_missing_region <- round(missing_region/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(!is.na(region)) 

# missing IMD & RUC
missing_arealevel <- sum(is.na(babycheck_analysis$imd))
percent_missing_arealevel <- round(missing_arealevel/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(!is.na(imd)) 

# missing first time mother 
missing_ftm <- sum(is.na(babycheck_analysis$first_time_mum))
percent_missing_ftm <- round(missing_ftm/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(!is.na(first_time_mum)) 

# Gypsy or Irish Traveller 
CM_goit <- sum(babycheck_analysis$M_eth18_2011 == "Gypsy or Irish Traveller" | babycheck_analysis$C_eth18_2011 == "Gypsy or Irish Traveller")
percent_CM_goit <- round(CM_goit/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(M_eth18_2011 != "Gypsy or Irish Traveller") %>% 
  filter(C_eth18_2011 != "Gypsy or Irish Traveller") 
  
# Arab
CM_arab <- sum(babycheck_analysis$M_eth18_2011 == "Arab" | babycheck_analysis$C_eth18_2011 == "Arab")
percent_CM_arab <- round(CM_arab/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(M_eth18_2011 != "Arab") %>% 
  filter(C_eth18_2011 != "Arab") 

# Born financial year 2005-06
yob2005 <- sum(babycheck_analysis$financialyear == "2005-06")
percent_yob2005 <- round(yob2005/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(financialyear != "2005-06") 

# follow up end before follow up start (dob)
fuend_before_fustart <- sum(babycheck_analysis$C_pc_fu_end < babycheck_analysis$dob)
percent_fuend_before_fustart <- round(fuend_before_fustart/total_child*100, 3)
sourcepop <- sourcepop %>% 
  filter(C_pc_fu_end >= dob)

# follow-up end before the end of 12 weeks (already counted in the section above)
babycheck_sensanalysis <- sourcepop %>% 
  filter(C_pc_fu_end >= dob %m+% weeks(13)) 

# not registered by 1 year
not_reg_by_1y <- sum(babycheck_analysis$C_regstartdate > babycheck_analysis$dob %m+% years(1)) 
percent_not_reg_by_1y <- round(not_reg_by_1y/total_child*100, 3)
babycheck_sensanalysis <- babycheck_sensanalysis %>% 
  filter(C_regstartdate <= dob %m+% years(1)) 

# not registered by 4 weeks 
not_reg_by_4wk <- sum(babycheck_analysis$C_regstartdate > babycheck_analysis$dob %m+% weeks(4)) 
percent_not_reg_by_4wk <- round(not_reg_by_4wk/total_child*100, 3)
babycheck_mainanalysis <- babycheck_sensanalysis %>% 
  filter(C_regstartdate <= dob %m+% weeks(4))

# final count of children 
n_sourcepop <- n_distinct(sourcepop$C_patid) # 1339494
n_babycheck_mainanalysis <- n_distinct(babycheck_mainanalysis$C_patid) 
n_babycheck_sensanalysis<- n_distinct(babycheck_sensanalysis$C_patid) 

# remove outcome data from source pop file and check no NAs
sourcepop <- sourcepop %>% 
  select(-c(babycheck, certainty, hp, f2f, obsdate, 
            maternalcheck, maternalcheck_4to12reg, maternalcheck_12m, 
            babycheck_maternalcheck, babycheck_maternalcheck_4to12reg, babycheck_maternalcheck_12m,
            babycheck_vax))
summary(is.na(sourcepop))

# save
save(sourcepop, file="~/sourcepop.Rdata")
save(babycheck_mainanalysis, file="~/babycheck_mainanalysis.Rdata")
save(babycheck_sensanalysis, file="~/babycheck_sensanalysis.Rdata")

# save counts
indeterminate_sex <- as.data.frame(indeterminate_sex)
percent_indeterminate_sex <- as.data.frame(percent_indeterminate_sex)
multiple_birth <- as.data.frame(multiple_birth)
percent_multiple_birth <- as.data.frame(percent_multiple_birth)

missing_region <- as.data.frame(missing_region)
percent_missing_region <- as.data.frame(percent_missing_region)
missing_arealevel <- as.data.frame(missing_arealevel)
percent_missing_arealevel <- as.data.frame(percent_missing_arealevel)
missing_ftm <- as.data.frame(missing_ftm)
percent_missing_ftm <- as.data.frame(percent_missing_ftm)
CM_goit <- as.data.frame(CM_goit)
percent_CM_goit <- as.data.frame(percent_CM_goit)
CM_arab <- as.data.frame(CM_arab)
percent_CM_arab <- as.data.frame(percent_CM_arab)
yob2005 <- as.data.frame(yob2005)
percent_yob2005 <- as.data.frame(percent_yob2005)

fuend_before_fustart <- as.data.frame(fuend_before_fustart)
percent_fuend_before_fustart <- as.data.frame(percent_fuend_before_fustart)
fuend_before12wk <- as.data.frame(fuend_before12wk)
percent_fuend_before12wk <- as.data.frame(percent_fuend_before12wk)
not_reg_by_1y <- as.data.frame(not_reg_by_1y)
percent_not_reg_by_1y <- as.data.frame(percent_not_reg_by_1y)
not_reg_by_4wk <- as.data.frame(not_reg_by_4wk)
percent_not_reg_by_4wk <- as.data.frame(percent_not_reg_by_4wk)

n_sourcepop <- as.data.frame(n_sourcepop)
n_babycheck_mainanalysis <- as.data.frame(n_babycheck_mainanalysis)
n_babycheck_sensanalysis <- as.data.frame(n_babycheck_sensanalysis)

babycheck_exclusions <- cbind(indeterminate_sex, percent_indeterminate_sex,
                              multiple_birth, percent_multiple_birth,
                              missing_region, percent_missing_region,
                              missing_arealevel, percent_missing_arealevel,
                              missing_ftm, percent_missing_ftm,
                              CM_goit, percent_CM_goit,
                              CM_arab, percent_CM_arab,
                              yob2005, percent_yob2005,
                              fuend_before_fustart, percent_fuend_before_fustart,
                              fuend_before12wk, percent_fuend_before12wk,
                              not_reg_by_1y, percent_not_reg_by_1y,
                              not_reg_by_4wk, percent_not_reg_by_4wk,
                              n_sourcepop, n_babycheck_mainanalysis, n_babycheck_sensanalysis)

write_csv(babycheck_exclusions, file = "~/babycheck_exclusions.csv")

# Create GP registration variable -----------------------------------------

load(file="~/babycheck_sensanalysis.Rdata")

# create age in weeks since GP registration variable
babycheck_sensanalysis <- babycheck_sensanalysis %>% 
  mutate(wk_since_reg = interval(dob, C_regstartdate)/weeks(1)) %>% 
  mutate(wk_since_reg = round(wk_since_reg, 0)) %>% 
  mutate(wk_since_reg = ifelse(wk_since_reg < 0, 0, wk_since_reg)) # for those registering before estimated birth date

# check ranges
range(babycheck_sensanalysis$wk_since_reg)
hist(babycheck_sensanalysis$wk_since_reg)

# save
save(babycheck_sensanalysis, file="~/babycheck_sensanalysis.Rdata")
