# MODEL FITTING: AGE DEPENDENT CHANGES IN INTRINSIC TIMESCALE -- ACF SUM

library(stringr)
library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(ggseg)
library(ggsegSchaefer)
library(paletteer)
library(pals)
library(ggseg3d)
library(ggplot2)
library(scales)
library(tibble)
library(purrr)

source("/Volumes/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/code/fcn_GAM_timescale.R")

# load and prepare data for GAMs
dataset <- 'HCPD_rest_age' # 'HBN_rest_noSI_age' 'HCPD_rest_age' 'HCPYA_rest'
metric <- 'timescale'
dsmethod <- 'TRorig'

if (metric == 'timescale'){
  metricFileName <- 'timescale_acf'}else{
    metricFileName <- metric}

if (metric == 'alff'){regionmeasure <- 'meanALFF'}
if (metric == 'falff'){regionmeasure <- 'meanfALFF'}
if (metric == 'timescale'){regionmeasure <- 'meanTS'}
if (metric == 'rad'){regionmeasure <- 'meanRAD'}
if (metric == 'hurst'){regionmeasure <- 'meanHurst'}
if (metric == 'ac1'){regionmeasure <- 'meanAC1'}

project_path <- '/Volumes/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/'
data_path <- paste(project_path, sprintf('data/%s/', metric), sep = "")
outpath <- paste(project_path, sprintf('results/%s/', metric), sep = "")
ts.schaefer400.all <- read.csv(paste(data_path, 
                                     sprintf('%s_acf/%s/concat/%s_concat_%s_%s_Schaefer_400-7_forR.tsv', 
                                             dataset, dsmethod, dataset, 
                                             metricFileName, dsmethod), 
                                     sep = ""), 
                               sep = '\t')
# will use dataframe as a covariate
ts.schaefer400.all$dataset <- as.factor(ts.schaefer400.all$study)
# will use sex as a covariate
ts.schaefer400.all$sex <- as.factor(ts.schaefer400.all$sex)

#################################
# plot Mean TS across individuals
#################################
corticalMean <- colMeans(ts.schaefer400.all[, 1:400], na.rm = TRUE)
corticalMean <- as.data.frame(corticalMean)
corticalMean <- tibble::rownames_to_column(corticalMean, "region")

# SA axis
sa.schaefer400 <- read.csv(paste(project_path, 
                                 'data/SchaeferParcellation/SArank_schaefer400_7Networks.csv', 
                                 sep = ""))
# sa.schaefer400 <- sa.schaefer400 %>% select(-X)
colnames(sa.schaefer400) <- c("SA.rank", "region")

corticalMean$region <- gsub("X", "", corticalMean$region)

corticalMean <- merge(corticalMean, sa.schaefer400, by = "region")

ggseg(.data = corticalMean, atlas = "schaefer7_400", 
      mapping=aes(fill = corticalMean, colour=I("#e9ecef"),
                  size=I(.03)), position = c("stacked")) + 
  theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::ocean.matter",  # grDevices::Grays
                                    na.value="transparent", 
                                    direction = -1, 
                                    limits = c(min(corticalMean$corticalMean),
                                               max(corticalMean$corticalMean)), 
                                    oob = squish)

ggsave(filename = paste(outpath, sprintf('%s_brainmap_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_brainmap_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)


# plt mean TS ranks on the brain
# Rank meanTS values
corticalMean$corticalMeanRank <- rank(
  corticalMean$corticalMean, 
  ties.method = "average"
)

ggseg(.data = corticalMean, atlas = "schaefer7_400", 
      mapping=aes(fill = corticalMeanRank, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", 
                                    na.value="transparent", 
                                    direction = -1, 
                                    limits = c(min(corticalMean$corticalMeanRank),
                                               max(corticalMean$corticalMeanRank)), 
                                    oob = squish)

ggsave(filename = paste(outpath, sprintf('%s_brainmap_%s_%s_ranks.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_brainmap_%s_%s_ranks.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

#################################
# Plot mean TS vs SA rank
#################################
# SA rank comparison
rho = cor.test(corticalMean$corticalMean, 
               corticalMean$SA.rank, method = c("spearman"))

ggplot(corticalMean, aes(x = SA.rank, y = corticalMean, 
                             fill = SA.rank)) + 
  geom_point(aes(color = SA.rank), shape = 21, size = 2) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
                       guide = "colourbar", aesthetics = "fill", name = NULL, 
                       midpoint = 200) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
                       guide = "colourbar", aesthetics = "color", name = NULL, 
                       midpoint = 200) +
  labs(x="\nS-A rank", y="mean TS\n") +
  ggtitle(rho$estimate) +
  geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("gray70"),.7), 
              col = "black", linewidth = 1) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12, family = "Arial", 
                                 color = c("black")), 
        axis.title = element_text(size = 12, family = "Arial", 
                                  color = c("black")))

ggsave(filename = paste(outpath, sprintf('%s_SArank_meanTS_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())
ggsave(filename = paste(outpath, sprintf('%s_SArank_meanTS_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())

#################################
# Plot mean TS vs SA rank for age bins
#################################
# --- Prepare SA rank table ---
sa.schaefer400 <- read.csv(paste(project_path, 
                                 'data/SchaeferParcellation/SArank_schaefer400_7Networks.csv', 
                                 sep = ""))
colnames(sa.schaefer400) <- c("SA.rank", "region")
sa.schaefer400 <- sa.schaefer400 %>%
  mutate(region = as.character(region))

# --- Helper: given a subject subset, compute Spearman rho across regions ---
compute_rho_for_group <- function(subdf) {
  # mean across subjects for each region (columns 1:400)
  cm <- colMeans(subdf[, 1:400, drop = FALSE], na.rm = TRUE)
  df_cm <- tibble(region = names(cm), corticalMean = as.numeric(cm)) %>%
    mutate(region = gsub("^X", "", region)) %>%
    left_join(sa.schaefer400, by = "region")
  
  # Spearman correlation across regions
  ct <- suppressWarnings(
    cor.test(df_cm$corticalMean, df_cm$SA.rank,
             method = "spearman", exact = FALSE)
  )
  
  tibble(
    rho = unname(ct$estimate),
    p   = ct$p.value,
    n_regions = sum(complete.cases(df_cm$corticalMean, df_cm$SA.rank))
  )
}

# --- Create quarter-year bins ---
breaks_q <- seq(8, 22, by = 0.25)

res_by_bin <- ts.schaefer400.all %>%
  filter(!is.na(age), age >= 8, age < 22) %>%
  mutate(
    # get integer bin index directly
    bin_id = cut(age, breaks = breaks_q, include.lowest = TRUE,
                 right = FALSE, labels = FALSE),
    # midpoint of the bin
    age_center = breaks_q[bin_id] + 0.125
  ) %>%
  group_by(age_center) %>%
  group_modify(~ compute_rho_for_group(.x) %>% mutate(n_subjects = nrow(.x))) %>%
  ungroup() %>%
  arrange(age_center)

# Inspect results
print(res_by_bin)
# A tibble with columns: age_bin, rho, p, n_regions, n_subjects

# --- Plot rho vs. age bin ---
rhoOfRho = cor.test(res_by_bin$age_center, 
                    res_by_bin$rho, method = c("spearman"))

ggplot(res_by_bin, aes(x = age_center, y = rho)) +
  geom_point(color = "#F9BE85", size = 2) +
  theme_classic() +
  geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("gray70"),.7), 
              col = "black", size = .25, linewidth = 1) +
  labs(x = "Age (quarter-year bin center)",
       y = expression(Spearman~rho~"(corticalMean × SA rank)")) +
  ggtitle(rhoOfRho$estimate) +
  # scale_x_continuous(breaks = seq(8, 22, by = 1))
  theme(
    axis.text = element_text(size=12, family = "Arial", 
                             color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", 
                                color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", 
                                color = c("black"))) +
  scale_y_continuous(
    limits = c(-0.1, 0.5),  # HBN: c(-0.2, 0.4)
    expand = c(0, .01)
  ) +
  if (dataset == 'HCPYA_rest'){
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36, 38), 
                       expand = c(0,.45))
  }else{
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                       expand = c(0,.45))
  }

ggsave(filename = paste(outpath, sprintf('%s_quarterage_SA_TS_%s_scatter_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())
ggsave(filename = paste(outpath, sprintf('%s_quarterage_SA_TS_%s_scatter_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())

# --- Bin ages per year and compute rho per bin ---
# Bin definition: [k, k+1) gets labeled with center k + 0.5
res_by_bin <- ts.schaefer400.all %>%
  filter(!is.na(age), age >= 8, age < 22) %>%
  mutate(age_center = floor(age) + 0.5) %>%
  group_by(age_center) %>%
  # skip empty/small bins if you want; otherwise compute blindly
  group_modify(~ compute_rho_for_group(.x) %>% mutate(n_subjects = nrow(.x))) %>%
  ungroup()

# Inspect results
print(res_by_bin)
# A tibble with columns: age_center, rho, p, n_regions, n_subjects

# --- Plot rho vs. age bin ---
rhoOfRho = cor.test(res_by_bin$age_center, 
                    res_by_bin$rho, method = c("spearman"))

ggplot(res_by_bin, aes(x = age_center, y = rho)) +
  geom_point(color = "#F9BE85", size = 2) +
  theme_classic() +
  geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("gray70"),.7), 
              col = "black", size = .25, linewidth = 1) +
  labs(x = "Age (one-year bin center)", 
       y = expression(Spearman~rho~"(corticalMean × SA rank)")) +
  ggtitle(rhoOfRho$estimate) + 
  # scale_x_continuous(breaks = seq(8.5, 21.5, by = 1))
  theme(
    axis.text = element_text(size=12, family = "Arial", 
                             color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", 
                                color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", 
                                color = c("black"))) +
  # scale_y_continuous(
  #   limits = c(-0.1, 0.5),  # HBN: c(0, 0.4)
  #   expand = c(0, .01)
  # ) +
  if (dataset == 'HCPYA_rest'){
      scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36, 38), 
                         expand = c(0,.45),
                         limits = c(22, 38))
  }else{
      scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                         expand = c(0,.45),
                         limits = c(8, 22))
  }

ggsave(filename = paste(outpath, sprintf('%s_yearlyage_SA_TS_%s_scatter_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())
ggsave(filename = paste(outpath, sprintf('%s_yearlyage_SA_TS_%s_scatter_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())

#################################
# fit GAMs across individuals for mean TS 
#################################
# k = 3, fixed edf
ctx.predicted.ts <- gam.smooth.predict(measure = 'ts', atlas = 'schaefer400', 
                                       dataset = 'all', 
                                       region = regionmeasure, 
                                       smooth_var = 'age', 
                                       covariates = 'sex + meanFD', 
                                       knots = 3, set_fx = FALSE, 
                                       increments = 300)
# get predicted.smooth df from function output
preddata <- as.data.frame(ctx.predicted.ts[3])

ggplot(data = ts.schaefer400.all, aes(x = age, y = meanTS)) +
  geom_point(color = "#F9BE85", size = 2) + # color = "#115c25"
  geom_ribbon(data = preddata, aes(x = age, y = .fitted,  
                                   ymin = .lower_ci, ymax = .upper_ci), 
              alpha = .7, linetype = 0, fill="gray70") +
  geom_line(data = preddata, aes(x = age, y = .fitted), 
            color = "black", linewidth = 1) +
  labs(x="\nage", y=paste(regionmeasure, "\n", sep = "")) +
  theme_classic() +
  theme(
    axis.text = element_text(size=12, family = "Arial", 
                             color = c("black")),
    axis.title.x = element_text(size=12, family ="Arial", 
                                color = c("black")),
    axis.title.y = element_text(size=12, family ="Arial", 
                                color = c("black"))) +
  if (dataset == 'HCPYA_rest'){
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36, 38), 
                       expand = c(0,.45))
  }else{
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                       expand = c(0,.45))
  }

ggsave(filename = paste(outpath, sprintf('%s_age_%s_scatter_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())
ggsave(filename = paste(outpath, sprintf('%s_age_%s_scatter_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())

# same plot with site info
if (dataset == 'HBN_rest_noSI_age'){
  ggplot(data = ts.schaefer400.all, aes(x = age, y = meanTS)) +
    geom_point(aes(color = study_site), size = 2) + # color = "#115c25"
    geom_line(data = preddata, aes(x = age, y = .fitted)) +
    geom_ribbon(data = preddata, aes(x = age, y = .fitted,
                                     ymin = .lower_ci, ymax = .upper_ci),
                alpha = .7, linetype = 0,) +
    labs(x="\nage", y=paste(regionmeasure, "\n", sep = "")) +
    theme_classic() +
    theme(
      axis.text = element_text(size=12, family = "Arial", 
                               color = c("black")),
      axis.title.x = element_text(size=12, family ="Arial", 
                                  color = c("black")),
      axis.title.y = element_text(size=12, family ="Arial", 
                                  color = c("black"))) +
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                       expand = c(0,.45))
  
  ggsave(filename = paste(outpath, sprintf('%s_age_%s_scatter_site_%s.png', 
                                           dataset, metric, dsmethod), 
                          sep = ""), 
         plot = last_plot())
  ggsave(filename = paste(outpath, sprintf('%s_age_%s_scatter_site_%s.svg', 
                                           dataset, metric, dsmethod), 
                          sep = ""), 
         plot = last_plot())
}
#################################
# Region-wise GAMs
#################################
# Fit GAM Models for each region (Age Effects)
gam.age.schaefer <- matrix(data=NA, nrow=400, ncol=10) 

# list of regions to run gam.fit.smooth function on below
schaefer.regions <- names(ts.schaefer400.all[0:400]) %>% 
  as.data.frame() %>% set_names("region") 

# for each schaefer region, run gam smooth function
for(row in c(1:nrow(schaefer.regions))){
  region <- schaefer.regions$region[row] 
  GAM.RESULTS <- gam.fit.smooth(measure = "ts", atlas = "schaefer400", 
                                dataset = "all", 
                                region = region, smooth_var = "age", 
                                covariates = "sex + meanFD", 
                                knots = 3, set_fx = FALSE, 
                                stats_only = FALSE)
  # append results to output df
  gam.age.schaefer[row,] <- GAM.RESULTS}

gam.age.schaefer <- as.data.frame(gam.age.schaefer)
colnames(gam.age.schaefer) <- c("region", "GAM.age.Fvalue", "GAM.age.pvalue",
                                "GAM.age.partialR2", "Anova.age.pvalue",
                                "age.onsetchange", "age.peakchange",
                                "minage.decrease", "maxage.increase",
                                "age.maturation")
cols = c(2:10)    
gam.age.schaefer[,cols] = apply(gam.age.schaefer[,cols], 2, 
                                function(x) as.numeric(as.character(x)))
write.csv(gam.age.schaefer, 
          paste(outpath, 
                sprintf('csvFiles/%s_%s_age_statistics_%s.csv', 
                        dataset, metric, dsmethod), 
                sep=""),
          row.names = F, quote = F)
rm(gam.age.schaefer)
gc()

#################################
# Read region-wise GAM results from above and save with another format
#################################
gam.age.schaefer <- read.csv(paste(outpath, 
                                   sprintf('csvFiles/%s_%s_age_statistics_%s.csv', 
                                                    dataset, metric, dsmethod), 
                                   sep=""))
# SA axis
sa.schaefer400 <- read.csv(paste(project_path, 
                                 'data/SchaeferParcellation/SArank_schaefer400_7Networks.csv', 
                                 sep = ""))
# sa.schaefer400 <- sa.schaefer400 %>% select(-X)
colnames(sa.schaefer400) <- c("SA.rank", "region")

gam.age.schaefer$region <- gsub("X", "", gam.age.schaefer$region)

gam.age.schaefer <- merge(gam.age.schaefer, sa.schaefer400, by = "region")

csvR2 <- data.frame(gam.age.schaefer$region)
csvR2$partialR2 <- gam.age.schaefer$GAM.age.partialR2
outputPath <- paste(outpath, sprintf('csvFiles/%s_%s_age_r2_%s.csv', 
                                     dataset, metric, dsmethod), 
                    sep="")
write.csv(csvR2, outputPath, row.names=FALSE)

#################################
# Plot region-wise GAM results
#################################
# Effect size: histogram
ggplot(gam.age.schaefer, aes(x = GAM.age.partialR2)) + 
  geom_histogram(binwidth=.01, fill="darkcyan", color="#e9ecef", alpha=0.9) + 
  theme_bw()
ggsave(filename = paste(outpath, sprintf('%s_histogram_partialR2_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_histogram_partialR2_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# Effect size: brain
maxval <- max(abs(gam.age.schaefer$GAM.age.partialR2))

ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = GAM.age.partialR2, colour=I("#e9ecef"),
                  size=I(.03)), position = c("stacked")) + 
  theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::warmcool", 
                                    na.value="transparent", direction = -1, 
                                    limits = c(-maxval, maxval), 
                                    oob = squish)
# grDevices::RdPu
# grDevices::PinkYl
# 
# pals::ocean.matter
ggsave(filename = paste(outpath, sprintf('%s_brainmap_partialR2_%s_%s_v2.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_brainmap_partialR2_%s_%s_v2.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# Effect size: ranks
# Replace GAM.age.partialR2 with ranks robustly
# (1) Rank the partial R2 values
gam.age.schaefer$GAM.age.rankR2 <- rank(
  gam.age.schaefer$GAM.age.partialR2, 
  ties.method = "average"
)

# (2) Find the nearest negative value, if it exists
neg_vals <- gam.age.schaefer$GAM.age.partialR2[
  gam.age.schaefer$GAM.age.partialR2 < 0
]

if (length(neg_vals) > 0) {
  # (3) Center ranks around the largest negative value (i.e., nearest to 0)
  nearestNeg <- max(neg_vals)
  nearestNegIdx <- which(gam.age.schaefer$GAM.age.partialR2 == nearestNeg)
  nearestNegRank <- gam.age.schaefer$GAM.age.rankR2[nearestNegIdx]
  
  gam.age.schaefer$GAM.age.rankR2 <- 
    gam.age.schaefer$GAM.age.rankR2 - (nearestNegRank + 1)
  
  centered <- TRUE
} else {
  # No negatives — do not shift the ranks
  centered <- FALSE
}

# (4) Max absolute value (for plotting scale)
maxval <- max(abs(gam.age.schaefer$GAM.age.rankR2), na.rm = TRUE)

# brain ranks
ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = GAM.age.rankR2, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent", 
                                    direction = -1, 
                                    limits = c(-maxval, maxval), 
                                    oob = squish) 

ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_rank_partialR2_%s_%s_v2.png', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_rank_partialR2_%s_%s_v2.svg', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# Effect size: brain p-values from Anova
pvalues = gam.age.schaefer$Anova.age.pvalue
pvaluesfdrs<-p.adjust(pvalues, method="BH")

Anovasignumber = sum(pvaluesfdrs < 0.05, na.rm=TRUE)
pvaluesfdrs[pvaluesfdrs >= 0.05] <- NA
gam.age.schaefer$Anova.age.pvaluefdr <- pvaluesfdrs

ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = Anova.age.pvaluefdr, colour=I("#e9ecef"), 
                  size=I(.03)), 
      position = c("stacked")) + theme_void() + ggtitle(Anovasignumber) +
  paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent", 
                                    direction = -1, 
                                    limits = c(0, 0.05), 
                                    oob = squish) 

ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_pval_partialR2_%s_%s_ANOVA.png', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_pval_partialR2_%s_%s_ANOVA.svg', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# Effect size: significant ranks
pvalues = gam.age.schaefer$Anova.age.pvalue
pvaluesfdrs <- p.adjust(pvalues, method="BH")
rankR2sig <- gam.age.schaefer$GAM.age.rankR2
rankR2sig[(pvaluesfdrs >= 0.05)] <- NA
gam.age.schaefer$GAM.age.rankR2sig <- rankR2sig

maxval <- max(abs(gam.age.schaefer$GAM.age.rankR2sig), na.rm=T)

ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = GAM.age.rankR2sig, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::warmcool", na.value="transparent", 
                                    direction = -1, 
                                    limits = c(-maxval, maxval), 
                                    oob = squish) 

ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_ranksig_partialR2_%s_%s_ANOVA.png', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, 
                        sprintf('%s_brainmap_ranksig_partialR2_%s_%s_ANOVA.svg', 
                                dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# save age results with ranks and fdr corrected pvalues
outputPath <- paste(outpath, sprintf('csvFiles/%s_%s_age_statistics_full_%s.csv', 
                                     dataset, metric, dsmethod), 
                    sep="")
write.csv(gam.age.schaefer, outputPath, row.names=FALSE)

#################################
# Plot region-wise GAM results: SA rank comparison
#################################
# SA rank comparison
rho = cor.test(gam.age.schaefer$GAM.age.partialR2, 
               gam.age.schaefer$SA.rank, method = c("spearman"))

ggplot(gam.age.schaefer, aes(x = SA.rank, y = GAM.age.partialR2, 
                             fill = SA.rank)) + 
  geom_point(aes(color = SA.rank), shape = 21, size = 2) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
                       guide = "colourbar", aesthetics = "fill", name = NULL, 
                       midpoint = 200) +
  scale_fill_gradient2(low = "goldenrod1", mid = "white", high = "#6f1282", 
                       guide = "colourbar", aesthetics = "color", name = NULL, 
                       midpoint = 200) +
  labs(x="\nS-A rank", y="partial R2\n") +
  ggtitle(rho$estimate) +
  geom_smooth(method = 'lm', se = TRUE, fill = alpha(c("gray70"),.7), 
              col = "black", linewidth = 1) +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 12, family = "Arial", 
                                 color = c("black")), 
        axis.title = element_text(size = 12, family = "Arial", 
                                  color = c("black")))

ggsave(filename = paste(outpath, sprintf('%s_SArank_partialR2_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())
ggsave(filename = paste(outpath, sprintf('%s_SArank_partialR2_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       plot = last_plot())

#################################
# Region-wise GAM fitted value predictions
#################################
# number of ages to predict timescale at
np <- 300

gam.smooths.schaefer <- matrix(data=NA, ncol=7) 
colnames(gam.smooths.schaefer) <- c("age", ".fitted", ".se",
                                    ".lower_ci", ".upper_ci", 
                                    "index", "region")

gam.peaks.schaefer <- matrix(data=NA, ncol=2) 
colnames(gam.peaks.schaefer) <- c("region", "age.peak")

# for each schaefer region, run the gam.smooth.predict function
for(row in c(1:nrow(schaefer.regions))){
  region <- schaefer.regions$region[row] 
  GAM.SMOOTH <- gam.smooth.predict(measure = "ts", atlas = "schaefer400", 
                                   dataset = "all", 
                                   region = region, smooth_var = "age", 
                                   covariates = "sex + meanFD", 
                                   knots = 3, set_fx = FALSE, increments = np)
  # get predicted.smooth df from function output
  preddata <- as.data.frame(GAM.SMOOTH[3])
  # region index
  preddata$index <- rep(x=row, np)
  # label
  preddata$region <- rep(x=GAM.SMOOTH[1], np)
  gam.smooths.schaefer <- rbind(gam.smooths.schaefer, preddata)
  # label and predicted timescale peak
  datapeak <- as.data.frame(cbind(GAM.SMOOTH[1], GAM.SMOOTH[2]))
  colnames(datapeak) <- c("region", "age.peak")
  gam.peaks.schaefer <- rbind(gam.peaks.schaefer, datapeak)
}

# remove empty initialization row
gam.smooths.schaefer <- gam.smooths.schaefer[-1,]
gam.smooths.schaefer$region <- as.character(gam.smooths.schaefer$region)
# remove empty initialization row
gam.peaks.schaefer <- gam.peaks.schaefer[-1,]
gam.peaks.schaefer$region <- as.character(gam.peaks.schaefer$region)
gam.peaks.schaefer$age.peak <- as.numeric(gam.peaks.schaefer$age.peak)
write.csv(gam.smooths.schaefer, 
          paste(outpath, 
                sprintf('csvFiles/%s_%s_age_predictedfits_%s.csv', 
                        dataset, metric, dsmethod), sep=""),
          row.names = F, quote = F) 
write.csv(gam.peaks.schaefer, 
          paste(outpath, sprintf('csvFiles/%s_%s_age_peaks_%s.csv', 
                                 dataset, metric, dsmethod), sep=""),
          row.names = F, quote = F)
rm(gam.smooths.schaefer)
gc()

#################################
# Region-wise GAM fitted value predictions: plot on brain
#################################
# reload results from above
gam.fitted.schaefer <- read.csv(paste(outpath, 
                                      sprintf('csvFiles/%s_%s_age_predictedfits_%s.csv', 
                                              dataset, metric, dsmethod), 
                                      sep="")) 
gam.fitted.schaefer$region <- gsub("X", "", gam.fitted.schaefer$region)
gam.fitted.schaefer <- inner_join(gam.fitted.schaefer, gam.age.schaefer, 
                                  by="region")

# brain min age onset change
ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = age.onsetchange, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + 
  theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", 
                                    na.value="transparent", 
                                    direction = 1, 
                                    limits = c(min(ts.schaefer400.all$age),
                                               max(ts.schaefer400.all$age)), 
                                    oob = squish) 

ggsave(filename = paste(outpath, sprintf('%s_ageonsetchange_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_ageonsetchange_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# brain peak change age
ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = age.peakchange, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + 
  theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", 
                                    na.value="transparent", 
                                    direction = 1, 
                                    limits = c(min(ts.schaefer400.all$age),
                                               max(ts.schaefer400.all$age)), 
                                    oob = squish) 

ggsave(filename = paste(outpath, sprintf('%s_agepeakchange_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_agepeakchange_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

# brain maturation age
ggseg(.data = gam.age.schaefer, atlas = "schaefer7_400", 
      mapping=aes(fill = age.maturation, colour=I("#e9ecef"), size=I(.03)), 
      position = c("stacked")) + 
  theme_void() + 
  paletteer::scale_fill_paletteer_c("pals::ocean.matter", 
                                    na.value="transparent", 
                                    direction = 1, 
                                    limits = c(min(ts.schaefer400.all$age),
                                               max(ts.schaefer400.all$age)), 
                                    oob = squish) 

ggsave(filename = paste(outpath, sprintf('%s_agematuration_%s_%s.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_agematuration_%s_%s.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)

#################################
# Region-wise smooth estimates: plot curves
#################################
np <- 300
gam.estimated.smooths.schaefer <- matrix(data=NA, ncol=4) 
colnames(gam.estimated.smooths.schaefer) <- 
  c("age",".estimate","index","region")

for(row in c(1:nrow(schaefer.regions))){
  region <- schaefer.regions$region[row] 
  GAM.ESTIMATES <- gam.estimate.smooth(measure = "ts", atlas = "schaefer400", 
                                       dataset = "all", 
                                       region = region, smooth_var = "age", 
                                       covariates = "sex + meanFD", 
                                       knots = 3, set_fx = FALSE, 
                                       increments = np)
  # region index
  GAM.ESTIMATES$index <- rep(x=row, np)
  # label
  GAM.ESTIMATES$region <- rep(x=region, np)
  gam.estimated.smooths.schaefer <- rbind(gam.estimated.smooths.schaefer, 
                                          GAM.ESTIMATES)
}

# remove empty initialization row
gam.estimated.smooths.schaefer <- gam.estimated.smooths.schaefer[-1,]

maxval <- max(abs(gam.age.schaefer$GAM.age.partialR2))

ggplot(gam.estimated.smooths.schaefer,
       aes(age, .estimate, group=index,
           color = gam.fitted.schaefer$GAM.age.partialR2)) +
  geom_line(size=.5, alpha = .8) + 
  paletteer::scale_color_paletteer_c("pals::warmcool", direction = -1,
                                     limits = c(-maxval, maxval),
                                     oob = squish) +
  theme_classic() +
  labs(x = "\nage", y = sprintf("%s (demeaned)\n", metric)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=12, family = "Arial", 
                                 color = c("black")), 
        axis.title = element_text(size=12, family = "Arial", 
                                  color = c("black"))) +
  if (dataset == 'HCPYA_rest'){
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36, 38), 
                       expand = c(0,.45))
  }else{
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                       expand = c(0,.45))
  }

ggsave(filename = paste(outpath, sprintf('%s_agefits_curves_%s_%s_v2.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_agefits_curves_%s_%s_v2.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)


# color model fits by SA rank
ggplot(gam.estimated.smooths.schaefer,
       aes(age, .estimate, group=index,
           color = gam.fitted.schaefer$SA.rank)) +
  geom_line(size=.5, alpha = .8) + 
  scale_color_gradient2(low = "goldenrod1",mid = "white", high = "#6f1282",
                        midpoint = 200, oob = scales::squish,
                        name = NULL, guide = "colourbar") +
  theme_classic() +
  labs(x = "\nage", y = sprintf("%s (demeaned)\n", metric)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=12, family = "Arial", 
                                 color = c("black")), 
        axis.title = element_text(size=12, family = "Arial", 
                                  color = c("black"))) +
  if (dataset == 'HCPYA_rest'){
    scale_x_continuous(breaks=c(22, 24, 26, 28, 30, 32, 34, 36, 38), 
                       expand = c(0,.45))
  }else{
    scale_x_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24), 
                       expand = c(0,.45))
  }

ggsave(filename = paste(outpath, sprintf('%s_agefits_curves_%s_%s_SArankcolorbar.png', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
ggsave(filename = paste(outpath, sprintf('%s_agefits_curves_%s_%s_SArankcolorbar.svg', 
                                         dataset, metric, dsmethod), 
                        sep = ""), 
       dpi = 300, width = 3 , height = 2)
