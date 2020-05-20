library(ggplot2)
source("scripts/plot_psmc.R")
source("scripts/theme_emily.R")
library(data.table)
library(plyr)
library(tidyr)
library(wesanderson)
library(scales)
options(scipen=999)
library(purrr)
source("scripts/theme_emily.R")
#install_github('tavareshugo/windowscanr')
library(dplyr)
library(windowscanr)

# Script to read in output from Macs simulation and plot results against sliding window heterozygosity
# Output from 4.3_macs_28_<sampleID>.sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#        Macs out         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Read in files

sim_files <- paste("data/macs/", list.files(path = "data/macs", pattern="*.txt"), sep = "")

# Read Macs output and extract number of segregating sites

read_sim <- function(x) {
  x <- readLines(x)
  seg_sites <- grep("^TOTAL_SITES", x) - 4
  seg_sites <- seg_sites
}

sims <- map(sim_files, read_sim)
macs_names <-  gsub(".txt", "", sim_files)
macs_names <-  gsub("data/macs/macs_", "", macs_names)
names(sims) <- macs_names

len <- 2500000

sims_df <- data.frame(V1 = unlist(sims)) %>%
  mutate(sim = rownames(.),
         het = V1 / len,
         id = "notrim") %>%
  select(sim, het, id) %>%
  separate(sim, c("trim", "id", "sim", "Sample"), sep = "_") %>%
  mutate(trim = as.numeric(trim),
         sim = as.numeric(sim)) %>%
  mutate(Sample_ID = case_when(Sample == "G449" ~ "G449",
                         Sample == "G445" ~ "G445",
                         Sample == "MSH009" ~ "#35552",
                         Sample == "MSH005" ~ "#34412",
                         Sample == "MSH638" ~ "#33556",
                         Sample == "MSH645" ~ "#36948"))
  

sims_df_ci <- sims_df %>%
  group_by(Sample_ID, trim) %>% 
  dplyr::summarise(CI_low = quantile(het, probs = c(0.025,0.975))[1],
                   CI_high = quantile(het, probs = c(0.025,0.975))[2])

sims_df <- left_join(sims_df, sims_df_ci, by = c("Sample_ID", "trim"))

mean_het <- sims_df %>%
  group_by(Sample_ID) %>%
  summarise(mean_het = mean(het))

ggplot(sims_df, aes(het)) +
  geom_histogram(aes(fill = trim)) + 
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, y = 20)) + 
  theme_emily() +
  geom_vline(aes(xintercept = mean_het), data = mean_het) + 
  facet_wrap(~Sample_ID, nrow = 3)


full <- filter(sims_df, trim == 28) %>%
  ungroup() %>%
  mutate(id = "sim") %>%
  #filter(Sample == "MSH009") %>%
  select(Sample_ID, het, id)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#      Sliding window heterozygosity      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Sliding window heterozygosity

win_size <- 2500000

het <- function(x){
  sum(x) / win_size
}


#~~ Read in filtered SNP file

snp_file <- fread("data/snps_focal/ORYX_geno_focal_biallelic_gdepth_miss_himaf.traw")


snps <- snp_file %>%
  select(-c(SNP, `(C)M`, COUNTED, ALT)) %>%
  pivot_longer(-c(CHR, POS), names_to = "Sample",
               values_to = "Geno") %>%
  mutate(Sample = gsub("_.+", "", Sample)) %>%
  filter(Geno == 1) %>%
  mutate(CHR = as.factor(CHR)) %>%
  mutate(Sample_ID = case_when(Sample == "G449" ~ "G449",
                             Sample == "G445" ~ "G445",
                             Sample == "MSH009" ~ "#35552",
                             Sample == "MSH005" ~ "#34412",
                             Sample == "MSH638" ~ "#33556",
                             Sample == "MSH645" ~ "#36948"))

# Run windowscanR
emp_het <- winScan(snps, position = "POS",
             values = "Geno",
             group = c("Sample_ID", "CHR"),
             win_size = 2500000, # 25 Mb
             win_step = 2500000,
             #funs = c("sum_na_rm", "len_na_rm")) %>%
             funs = "het")

#~~
emp <- emp_het %>%
  select(Sample_ID, Geno_het) %>%
  mutate(id = "emp")

colnames(emp) <- c("Sample_ID", "het", "id")
mean(emp$het)

df <- rbind(emp, full) %>%
  mutate(id = factor(id, levels = c("emp", "sim")),
         Sample_ID = factor(Sample_ID, levels = c("#33556", "#36948",
                                                  "#34412", "#35552",
                                                  "G445", "G449")))

df_ci <- filter(df, id == "sim") %>%
  dplyr::group_by(Sample_ID, id) %>% 
  dplyr::summarise(CI_low = quantile(het, probs = c(0.025,0.975))[1],
                   CI_high = quantile(het, probs = c(0.025,0.975))[2])

df_ci$Sample_ID <- factor(df_ci$Sample_ID, levels = c("#33556", "#36948",
                                                      "#34412", "#35552",
                                                      "G445", "G449"))

df_mean <- filter(df, id == "emp") %>%
  group_by(Sample_ID, id) %>% 
  dplyr::summarise(mean_het = mean(het))

df_mean$Sample_ID <- factor(df_mean$Sample_ID, levels = c("#33556", "#36948",
                                                      "#34412", "#35552",
                                                      "G445", "G449"))

#~~ Plotting

cbPalette <- wes_palette("Moonrise2")
colpal <- c("grey40", cbPalette[2])


png(file="figs/Macs_sim.png", units = "in", res = 300, height=9, width=8)

ggplot(df) +
  geom_histogram(aes(x = het, fill = id, col = id), alpha = 0.4, 
                 position = "identity", bins = 50) +
  geom_vline(data=df_mean, aes(xintercept=mean_het),
             lty = "dashed", col = "grey40",
             size = 0.7) +
  geom_errorbarh(data = df_ci,
                 aes(xmin = CI_low , xmax = CI_high, y = 230),
                 size = 0.7, color = cbPalette[2],
                 lty = "solid", height = 0) +
  facet_wrap(~Sample_ID, nrow = 3,
             ncol = 2) +
  scale_fill_manual(values = colpal,
                    name = "",
                    labels = c("Empirical", "Simulated")) +
  scale_colour_manual(values = colpal,
                      name = "",
                      labels = c("Empirical", "Simulated")) +
  theme_emily() +
  theme(legend.position=c(.8, .9)) +
  ylab("Count") + xlab("Genome-wide heterozygosity")

dev.off()


