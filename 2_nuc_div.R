library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(readxl)
library(plyr)
source("scripts/theme_emily.R")
library(wesanderson)
library(forcats)
library(ggrepel)
library(gridExtra)
options(scipen=999)
library(scales)
library(devtools)
devtools::install_github("thomasp85/patchwork")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Individual genome-wide heterozygosity  #  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Generated with one ind SFS using ANGSD for each chr
# run_5.1_angsdHet.sh &	run_5.2_realSFS.sh from https://github.com/elhumble/SHO_reseq_2020

import_files <- function(file){
  fread(file)
}

# Load ml estimates for each ind and each chr

ml_files <- paste("data/ind_het_ml/", list.files(path = "data/ind_het_ml", pattern="*.ml"), sep = "")
ml <- lapply(ml_files, import_files)

dataset_names <- list.files(path = "data/ind_het_ml", pattern = "*ml")
dataset_names <- gsub("_est.ml", "", dataset_names)
names(ml) <- dataset_names

# Summarise by individual

ind_het_chr <- ldply(ml, .id = "ID") %>%
  separate(ID, c("ID", "CHR")) %>%
  filter(CHR != 23) %>%
  dplyr::group_by(ID, CHR) %>%
  dplyr::summarise(V1 = sum(V1),
                   V2 = sum(V2)) %>%
  mutate(het = V2 / (V1 + V2)) %>%
  mutate(het = V2 / (V1 + V2),
         pop = case_when(ID = grepl("^G", ID) | grepl("^W", ID) ~ "EAD",
                         ID = grepl("^MSH2", ID) | grepl("^MSH3", ID) ~ "EAD",
                         ID = grepl("^MSH0", ID) ~ "EEP",
                         TRUE ~ "SSP")) %>%
  mutate(pop = as.factor(pop))


# 6 individuals

cbPalette <- wes_palette("Moonrise2")
cbPalette <- cbPalette[c(3,1,2)]

ind_het_6_chr <- filter(ind_het_chr, ID == "G449" | ID == "G445" | 
                      ID == "MSH009" | ID == "MSH005" | 
                      ID == "MSH638" | ID == "MSH645") %>%
  mutate(ID2 = case_when(ID == "G449" ~ "EAD 1",
                         ID == "G445" ~ "EAD 2",
                         ID == "MSH009" ~ "EEP 1",
                         ID == "MSH005" ~ "EEP 2",
                         ID == "MSH638" ~ "SSP 1",
                         ID == "MSH645" ~ "SSP 2")) %>%
  mutate(ID3 = case_when(ID == "G449" ~ "G449",
                         ID == "G445" ~ "G445",
                         ID == "MSH009" ~ "#35552",
                         ID == "MSH005" ~ "#34412",
                         ID == "MSH638" ~ "#33556",
                         ID == "MSH645" ~ "#36948"))
  

ind_het_6 <- ind_het_6_chr %>%
  dplyr::group_by(ID2, ID, pop) %>%
  dplyr::summarise(V1 = sum(V1),
            V2 = sum(V2),
            het = V2 / (V1 + V2))

summary(ind_het_6$het)

png(file="figs/genome_wide_het_pop_6.png", units = "in", res = 300, height=6, width=7)

ggplot(ind_het_6, aes(fct_reorder(pop, het), het, fill = pop)) +
  geom_boxplot() +
  theme_emily() +
  scale_fill_manual(values = cbPalette, name = "Population") +
  ylab("Genome-wide heterozygosity") + xlab("Population")

dev.off()

ind_het_6_chr <- ind_het_6_chr %>%
  mutate(ID3 = factor(ID3, levels = c("G449", "G445",
                                    "#35552", "#34412",
                                    "#33556", "#36948")))


het_6 <- ggplot(ind_het_6_chr, aes(ID3, het, fill = pop)) +
  geom_boxplot(outlier.size = -1, col = "grey40", alpha = 0.9) +
  geom_jitter(size = 2, shape = 21, alpha = 0.3, width = 0.1) +
  theme_emily() +
  scale_fill_manual(values = cbPalette[c(2,3,1)], name = "Population") +
  ylab("Genome-wide heterozygosity") + xlab("Individual ID") +
  theme(legend.position=c(.9, .2)) +
  ggtitle("B")
# +
#  geom_signif(comparisons = list(c("G445", "MSH009")), 
#              map_signif_level=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   Across species comparison   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ~ 10,000 oryx on the planet
# ~ 8,000 in texas ranches
# ~ 2-300 in Europe, AZA captive breeding
# ~ couple of thousand in UAE

# Read in species hets from Ekblom paper

sp_het <- read_excel("data/additional_sp_het.xlsx") %>%
  mutate(Common_Name = gsub("Bottlenose Dolphin", "Dolphin", Common_Name))

# Add oryx data based on 6 focal individuals

oryx <- data.frame(Species = "Oryx dammah",
                   Common_Name = "Scimitar-horned oryx",
                   Observed_Het = mean(ind_het_6$het),
                   Min_Het = min(ind_het_6$het),
                   Max_Het = max(ind_het_6$het),
                   N_census = 15000,
                   IUCN = "EW",
                   IUCN_grp = "Scimitar-horned oryx",
                   Label = "Scimitar-horned oryx",
                   Reference = "This_paper")

sp_het <- rbind(sp_het, oryx)

# Heterozygosity by IUCN status

ggplot(sp_het, aes(IUCN, Observed_Het)) +
  geom_boxplot() +
  theme_emily()

ggplot(sp_het, aes(IUCN_grp, Observed_Het)) +
  geom_boxplot() +
  theme_emily()

# Heterozygosity with census size

options(scipen=999)

# Plot

sp_het <- sp_het %>%
  mutate(IUCN = case_when(IUCN == "CR" ~ "CR",
                   TRUE ~ IUCN)) %>%
  mutate(IUCN = factor(IUCN, levels = c("EW", "CR", "EN", 
                                        "VU", "NT", "LC")))

print(levels(sp_het$IUCN))

sp_het_plot <- ggplot(sp_het, aes(N_census, Observed_Het, 
                                  fill = IUCN, label = Common_Name)) +
  geom_text_repel(data = subset(sp_het, N_census == 600000), # Dolphin
                  nudge_x = 0.5,
                  nudge_y = 0.00015, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 150000), # W Gorilla
                  nudge_x = 0.6,
                  nudge_y = -0.0001, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 90000), # Chimp
                  nudge_x = 0.7, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 50000), # Wolf
                  nudge_x = -0.5, 
                  nudge_y = 0.0002, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 30000), # Lion
                  nudge_x = 0.2, 
                  nudge_y = -0.0002, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, Common_Name == "Bonobo"), # Bonobo
                  nudge_x = -0.35,
                  nudge_y = -0.0001, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, Common_Name == "Scimitar-horned oryx"), # SHO
                  nudge_x = -0.9, segment.size = 0.2, fontface = "bold") +
  
  geom_text_repel(data = subset(sp_het, Common_Name == "Sumatran Orangutan"), # Sumatran Orang
                  nudge_x = -0.5,
                  nudge_y = 0.0003, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 8000), # Lynx
                  nudge_x = 0.5,
                  nudge_y = 0.0001, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 6700), # Cheetah
                  nudge_x = 0.1,
                  nudge_y = -0.0001, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census == 850), # Wolverine
                  nudge_x = -0.5,
                  nudge_y = 0.00015, segment.size = 0.2) +
  geom_text_repel(data = subset(sp_het, N_census < 850), # Baiji, Lynx, Tiger 
                  nudge_x = -0.6,
                  nudge_y = 0.0001, segment.size = 0.2) +
  geom_pointrange(size = 0.5, aes(ymin=Min_Het, ymax=Max_Het)) +
  geom_point(shape = 21, size = 3) +
  scale_fill_manual(values=c("black", "#CB2314", "#E58601", "#F2AD00", "#9C964A", "#0B775E"),
                    limits=c("EW", "CR", "EN",
                             "VU", "NT", "LC")) +
  theme_emily() +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1)) + 
  scale_x_continuous(trans='log10', 
                     breaks = c(10, 1000, 100000),
                     label = comma) +
  ylab("Genome-wide heterozygosity") + xlab("Census population size") +
  ggtitle("A")

dev.off()


# Arrange plots

png(file="figs/Figure_2.png", units = "in", res = 300, height=5, width=11)

grid.arrange(sp_het_plot, het_6,
             ncol = 2, nrow = 1,
             widths = c(2.4, 1.6))
dev.off()


#~~~~~~~~~~~~~~~~#
#  Long-term Ne  #
#~~~~~~~~~~~~~~~~#

mean(ind_het_6$het)
mean(ind_het_6$het) / (4 * 1.1e-08)

