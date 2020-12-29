#IL-6 inhibition and renal effects: A Mendelian Randomisation Study 
#Pre-publication reference: https://www.authorea.com/users/365323/articles/485466-inhibition-of-interleukin-6-signalling-and-renal-function-a-mendelian-randomization-study

#Import libraries 
library(dplyr)
#install.packages("forestplot")
library("forestplot")
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#install.packages("MendelianRandomization")
library(MendelianRandomization)

##Specify the working directory.
setwd("/Users/davidryan/Library/Mobile Documents/com~apple~CloudDocs/Documents/github/il6_mendelian_randomisation")
getwd()

#Exposure Data from UK Biobank 
#Previously selected as significant SNP (p<5x10-8) independently (clumped at 0.01) associated with CRP levels 

il6 <- read_exposure_data(
  filename = "ukbb_iv.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  eaf_col = "EAF", 
  pval_col = 'p-val')

#Outcome Data from CKDGen Consortium (Wuttke 2019 et al) - available from CKDGen website 

egfr <- read_outcome_data(
  filename = "egfr.csv",
  sep = ",",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2", 
  eaf_col = "Freq1",
  pval_col = "P-value")

bun <- read_outcome_data(
  filename = "bun.csv",
  sep = ",",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2", 
  eaf_col = "Freq1",
  pval_col = "P-value")

ckd <- read_outcome_data(
  filename = "ckd.csv",
  sep = ",",
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2", 
  eaf_col = "Freq1",
  pval_col = "P-value")

#Harmonise data 
il6_egfr <- harmonise_data(il6, egfr, action = 2)
il6_ckd <- harmonise_data(il6, ckd, action=2)
il6_bun <- harmonise_data(il6, bun, action=2)

#Primary Analysis
il6_egfr_res <- mr(il6_egfr, method_list=c("mr_ivw_mre","mr_simple_median", "mr_weighted_median", "mr_egger_regression"))
il6_bun_res <- mr(il6_bun, method_list=c("mr_ivw_mre","mr_simple_median", "mr_weighted_median", "mr_egger_regression"))
il6_ckd_res <- mr(il6_ckd, method_list=c("mr_ivw_mre","mr_simple_median", "mr_weighted_median", "mr_egger_regression"))

#Individual SNPs 
il6_egfr_single_snp <-  mr_singlesnp(il6_egfr,single_method = "mr_wald_ratio",all_method = c("mr_ivw", "mr_egger_regression"))
il6_bun_single_snp <-  mr_singlesnp(il6_bun,single_method = "mr_wald_ratio",all_method = c("mr_ivw", "mr_egger_regression"))
il6_ckd_single_snp <-  mr_singlesnp(il6_ckd,single_method = "mr_wald_ratio",all_method = c("mr_ivw", "mr_egger_regression"))

#Individual Forest Plots
mr_forest_plot(il6_egfr_single_snp , exponentiate = FALSE)
mr_forest_plot(il6_bun_single_snp , exponentiate = FALSE)
mr_forest_plot(il6_ckd_single_snp , exponentiate = TRUE)


#Results - CKD
ckd_labels <- c("IVW", "Simple Median", "Weighted Median", "MR-egger")
ckd_OR <- c(0.948, 0.922, 0.943, 1.013)
ckd_lower <- c(0.822, 0.754, 0.784, 0.720)
ckd_upper <- c(1.094, 1.128, 1.133, 1.425)

#bun
bun_labels <- c("IVW", "Simple Median", "Weighted Median", "MR-egger")
bun_estimate <- c(0.009, 0.005, 0.012, 0.004)
bun_lower <- c(-0.003, -0.017, -0.006, -0.023)
bun_upper <- c(0.021, 0.027, 0.030, 0.031)

#egfr 
egfr_labels <- c("IVW", "Simple Median", "Weighted Median", "MR-egger")
egfr_estimate <- c(0.001, 0.002, 0.002, 0.001)
egfr_lower <- c(-0.004, -0.007, -0.005, -0.012)
egfr_upper <- c(0.007, 0.010, 0.009, 0.014)


tiff("ckd_forest.tiff")
# make plot
fplot <- forestplot(labeltext = ckd_labels,
                    boxsize=0.15,
                    mean = ckd_OR, 
                    lower = ckd_lower, 
                    upper = ckd_upper, 
                    hrzl_lines=TRUE,
                    xlab= "Odds Ratio", 
                    title="IL6 inhibition and risk of CKD", 
                    zero=1, 
                    xticks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75),
                    graphwidth = "auto", 
                    lineheight = "auto")
# save plot

dev.off()

tiff("bun_forest.tiff")
# make plot
fplot <- forestplot(labeltext = bun_labels,
                    boxsize=0.15,
                    mean = bun_estimate, 
                    lower = bun_lower, 
                    upper = bun_upper, 
                    hrzl_lines=TRUE,
                    xlab= "Estimate", 
                    title="IL6 inhibition and Blood Urea Nitrogen", 
                    zero=0, 
                    xticks = c(-0.075, -0.05, -0.025, 0, 0.025, 0.05, 0.075),
                    graphwidth = "auto", 
                    lineheight = "auto")

# save plot

dev.off()


tiff("egfr_forest.tiff")

# make plot
fplot <- forestplot(labeltext = egfr_labels,
                    boxsize=0.15,
                    mean = egfr_estimate, 
                    lower = egfr_lower, 
                    upper = egfr_upper, 
                    hrzl_lines=TRUE,
                    xlab= "Odds Ratio", 
                    title="IL6 inhibition and eGFR", 
                    zero=0, 
                    xticks = c(-0.05, -0.025, 0, 0.025, 0.05),
                    graphwidth = "auto", 
                    lineheight = "auto")

# save plot

dev.off()

#Heterogeneity Testing 
mr_heterogeneity(il6_egfr)
mr_heterogeneity(il6_bun)
mr_heterogeneity(il6_ckd)








