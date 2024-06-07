library(tidyverse)
library(here) 
library(rmarkdown)
library(knitr)

library(kableExtra)
library(janitor)
library(scales)
library(ggpubr)
library(Rsubread)

library(pheatmap)
library(edgeR)

# BAM to count file
aligned <- featureCounts(here::here("data", "HG00097.mapped.ILLUMINA.bwa.GBR.exome.20130415.bam"), annot.inbuilt = "hg19", isPairedEnd = TRUE)

# Load data to dataframes
eDat <- read.delim(here::here("data", "GSE157103_formatted_eDat.txt"), sep = "\t")
pDat <- read.delim(here::here("data", "GSE157103_formatted_pDat.txt"), sep = "\t")

# Read dataframes
as_tibble(pDat)

# Transpose eDat dataframes so that eDat and pDat have same ID columns them match these IDs
eDat <- eDat %>%column_to_rownames(var = "gene")
all(colnames(eDat) == pDat$ID)

# Separation of patients by factors.
pDat %>%dplyr::count(COVID, Sex)
pDat %>%dplyr::count(COVID, ICU)
pDat <- pDat %>%mutate(Age_bracket = case_when(
                       Age < 20 ~ "below_20",
                       between(Age, 21, 40) ~ "21-40",
                       between(Age, 41, 60) ~ "41-60",
                       Age > 60 ~ "above_60"))
pDat$Age_bracket <- as.factor(pDat$Age_bracket)
pDat %>% count(COVID, Age_bracket)

# Plot graphs to visualize these separation
pDat %>% 
  ggplot(aes(x = Age_bracket, fill = Age_bracket)) +
  geom_bar(stat = "count") +
  facet_grid(~COVID)

# Plot this distribution of data by age bracket
pDat$Age_bracket <- fct_relevel(pDat$Age_bracket, c("below_20", "21-40", "41-60", "above_60"))  
pDat %>% 
  ggplot(aes(x = Age_bracket, fill = Age_bracket)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  facet_grid(~COVID) +
  labs(title = "Distribution of COVID Patients by Age Bracket")

# Make a new variable containing only protein datas and plot on a graph
proteins <- pDat %>% dplyr::select(COVID, Ferritin_ng.ml, CRP_mg.l, Procalcitonin_ng.ml, Lactate_mmol.l, Fibrinogen_mg.dL)
proteins <- proteins %>% pivot_longer(cols = 2:6, names_to = "protein", values_to = "measurement")
proteins %>% 
  ggplot(aes(x = COVID, y = measurement, fill = COVID)) +
  geom_boxplot() +
  scale_fill_manual(values = c("orange", "grey")) +
  facet_wrap(~protein, scales = "free") +
  theme_minimal()

# Plot heatmap to visualize data correlations by genders
pDat <- pDat %>% column_to_rownames(var = "ID")
h1 <- samp_cor %>% 
  pheatmap(clustering_distance_cols = "euclidean", clustering_method = "complete", cluster_rows = TRUE,
           show_colnames = FALSE, show_rownames = FALSE, 
           annotation_row = pDat[c("COVID", "Sex")], annotation_col = pDat[c("COVID", "Sex")], 
           main = "Sample Correlations")

# QC step: keeping only sequences with RPM >= 1 in all samples
e_fil <- eDat %>% 
  rownames_to_column(var = "gene") %>% 
  filter_at(vars(-gene), all_vars(. >= 1))  %>% 
  column_to_rownames(var = "gene")

# Normalize data by relative-Log Expression method
genes <- as.data.frame(row.names(e_fil))
norm <- DGEList(counts = e_fil2, samples = pDat, genes = genes, group = pDat$COVID)
eNormz <- calcNormFactors(norm, method = "RLE") 
eNormz <- cpm(eNorm)
eNormz <- as.data.frame(eNorm)

# Loading eNorm
eNorm <- read.delim(here::here("data", "eNorm.txt"), sep = "\t")
eNorm <- eNorm %>% column_to_rownames(var = "gene")
pDat <- read.delim(here::here("data", "GSE157103_formatted_pDat.txt"), sep = "\t")
pDat <- pDat %>% column_to_rownames(var = "ID")

# Transforming eNorm values to log2(x)+1 and transposing the result, then caculate PCA.
e_log2 <- log2(eNorm + 1)
t_log2 <- as.data.frame(t(e_log2))
pca <- prcomp(t_log2, scale = FALSE, center = TRUE)
summary <- data.frame(PC = 1:126, var_explained = (pca$sdev)^2 / sum((pca$sdev)^2), 
                      cumulative = cumsum(pca$sdev^2 / sum(pca$sdev^2)))
summary <- summary %>% mutate(cumulative_perc = cumulative*100)
summary <- summary[1:30,]
scores <- as.data.frame(pca$x)
scores <- scores[c(1:30)]

# Adding PCA scores into pDat dataFrames.
mDat <- cbind(pDat, scores)

# variation analysis 
variable_variance <- lmmatrix(dep = scores, ind = pDat[c(2:15)], metric = "Pvalue")
vv_plot <- variable_variance %>% as.data.frame() 
vv_plot <- as.data.frame(t(vv_plot))
vv_plot <- vv_plot %>% 
  mutate(Principle_Component = 1:30) %>% 
  dplyr::select(Principle_Component, everything())
vv_plot <- vv_plot %>% 
  pivot_longer(cols = -c(Principle_Component), names_to = "variables", values_to = "pval") 
vv_plot <- vv_plot %>% 
  mutate(pval_cat = case_when(
    pval > 0.05  ~ "> 0.05",
    pval < 0.05 & pval > 0.01 ~ "< 0.05",
    pval < 0.01 & pval > 0.001 ~ "< 0.01",
    pval < 0.001 ~ "< 0.001"))
vv_plot %>% 
  ggplot(aes(x = Principle_Component, y = variables, fill = pval_cat)) +
  geom_tile() + 
  theme_bw() +
  labs(x = "PC", y = "Variables" , fill = "P value")
