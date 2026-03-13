#**Meta-analysis: tACS session dose parameters - RAVLT**

# Date: March 10, 2025
# Author: Melody Marie Lawas
# Description: This R script performs a meta-analysis on the studies measuring
# RAVLT-IR change scores. Version 2. 
# Note: setwd("~/Practical_Project/Data Extraction")

library(readxl)
library(tidyverse)
library(metafor)   # escalc
library(meta)      # metagen, metacont, forest, metainf

library(dmetar)    # optional helper functions

# 1. Load Data ------------------------------------------------------------
file_path <- "PP_Outcome_Extraction.xlsx"
sheets <- readxl::excel_sheets(file_path)
df_raw <- readxl::read_excel(file_path, sheet = "MoCA")

# 2. Pre-processing & SD Calculation ---------------------------------------
# We need Mean_change and SD_change for EVERY study.
# For studies missing SD_change, we calculate it using r = 0.5

r_assumed <- 0.3 # change to check for sensitivity

df_clean <- df_raw %>%
  mutate(
    # Ensure numeric types
    across(starts_with(c("Mean", "SD", "Sample")), as.numeric),
    
    # Calculate Mean Change if only Pre/Post are available
    Mean_change_int = ifelse(is.na(Mean_change_int), (Mean_post_immediate_int - Mean_pre_int), Mean_change_int),
    Mean_change_con = ifelse(is.na(Mean_change_con), (Mean_post_immediate_con - Mean_pre_con), Mean_change_con),
    
    # Calculate SD Change using the formula: SD_change = sqrt(SD1^2 + SD2^2 - 2*r*SD1*SD2)
    SD_change_int = ifelse(is.na(SD_change_int), 
                           sqrt(SD_pre_int^2 + SD_post_immediate_int^2 - (2 * r_assumed * SD_pre_int * SD_post_immediate_int)), 
                           SD_change_int),
    SD_change_con = ifelse(is.na(SD_change_con), 
                           sqrt(SD_pre_con^2 + SD_post_immediate_con^2 - (2 * r_assumed * SD_pre_con * SD_post_immediate_con)), 
                           SD_change_con),
    
    study_label = paste0(Authors, " (", Year, ")")
  )

# 3. Calculate Effect Sizes (Hedges' g) ------------------------------------
# We use measure = "SMD" because we are comparing the CHANGE of Int vs CHANGE of Con
es_all <- escalc(
  measure = "SMD",
  m1i = Mean_change_int, sd1i = SD_change_int, n1i = Sample_Size_int,
  m2i = Mean_change_con, sd2i = SD_change_con, n2i = Sample_Size_con,
  data = df_clean,
  append = TRUE
)

# 4. Meta-Analysis ---------------------------------------------------------
meta_final <- metagen(
  TE = yi,
  seTE = sqrt(vi),
  studlab = study_label,
  data = es_all,
  sm = "SMD",
  common = FALSE,
  random = TRUE,
  method.tau = "REML",
  method.random.ci = "HK" # Hartung-Knapp
)

# Check if Study_Design (0 vs 1) actually impacts the effect size
meta_subgroup <- update(meta_final, subgroup = Study_Design)
print(meta_subgroup)

# 5. Forest Plot & Sensitivity -------------------------------------------------
summary(meta_final)

forest(meta_final,
       sortvar = yi,
       prediction = TRUE,
       xlab = "Standardized Mean Change Difference (Hedges' g)",
       leftcols = c("studlab", "effect", "ci", "w.random"),
       leftlabs = c("Study", "SMD", "95% CI","Weight"),
       rightcols = FALSE,
       colgap.forest.left = "0.5cm")

# Forest Plot with Subgroups
forest(meta_subgroup,
       sortvar = yi,
       prediction = TRUE, 
       print.subgroup.labels = TRUE,
       leftcols = c("studlab", "effect", "ci", "w.random"),
       leftlabs = c("Study", "Hedges' g", "95% CI", "Weight"),
       rightcols = FALSE,
       main = "RAVLT-IR Meta-Analysis by Study Design",
       col.subgroup = "blue",
       colgap.forest.left = "1cm")

# Standard Funnel Plot
funnel(meta_final, 
       xlab = "Hedges' g", 
       main = "Funnel Plot for RAVLT-IR")

# Leave-one-out to see if the "change-only" study is an outlier
inf <- metainf(meta_final)
forest(inf, xlab = "Leave-one-out Hedges' g")
