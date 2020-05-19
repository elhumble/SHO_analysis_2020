library(ggplot2)
source("scripts/plot_psmc.R")
source("scripts/theme_emily.R")
library(data.table)
library(plyr)
library(tidyr)
library(wesanderson)
library(scales)
options(scipen=999)

#~~ Specify variables for plotting

i.iteration=20 # Number of iterations to use in file
s=100 # bin size

#~~ Read in main PSMC files. Generated using 4.1_PSMC.sh & run_4.2_boot_PSMC.sh from https://github.com/elhumble/SHO_reseq_2020

psmc_files <- paste("data/psmc/", list.files(path = "data/psmc", pattern="*.psmc"), sep = "")

#~~ Run for multiple values of mu and g:

lo_mu <- 0.8e-8
mid_mu <- 1.1e-8
hi_mu <- 1.3e-8

lo_g <- 3
mid_g <- 6.2
hi_g <- 10

psmc_a <- lapply(psmc_files, psmc.result, i.iteration, mu = lo_mu, s, g = lo_g) # low
psmc_b <- lapply(psmc_files, psmc.result, i.iteration, mu = mid_mu, s, g = mid_g) # reasonable
psmc_c <- lapply(psmc_files, psmc.result, i.iteration, mu = hi_mu, s, g = hi_g) # high

#~~ Combine all results

dataset_names <-list.files(path = "data/psmc/", pattern="*.psmc")
dataset_names <-  gsub(".psmc", "", dataset_names)

names(psmc_a) <- dataset_names
names(psmc_b) <- dataset_names
names(psmc_c) <- dataset_names

psmc_a <- ldply(psmc_a, .id = "Sample")  
psmc_b <- ldply(psmc_b, .id = "Sample")  
psmc_c <- ldply(psmc_c, .id = "Sample")  


psmc_df <- rbind(psmc_a, psmc_b, psmc_c) %>%
  mutate(pop = ifelse(grepl("^G", Sample), "EAD",
                      ifelse(grepl("^MSH0", Sample), "EEP", "SSP"))) %>%
  mutate(Sample_ID = case_when(Sample == "G449" ~ "G449",
                               Sample == "G445" ~ "G445",
                               Sample == "MSH009" ~ "#35552",
                               Sample == "MSH005" ~ "#34412",
                               Sample == "MSH638" ~ "#33556",
                               Sample == "MSH645" ~ "#36948"))


#~~ Get bootstraps

boot_files <- paste("data/psmc_boot/", list.files(path = "data/psmc_boot", pattern="*boot*"), sep = "")
boot <- lapply(boot_files, psmc.result, i.iteration, mid_mu, s, mid_g)

dataset_names <-list.files(path = "data/psmc_boot/", pattern="*boot*")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(boot) <- dataset_names

boot_df <- ldply(boot, .id = "ID") %>%
  separate(ID, c("Sample", "Boot"), sep = "_boot", remove = F) %>%
  mutate(Boot = gsub("_", "", Boot)) %>%
  mutate(pop = ifelse(grepl("^G", Sample), "EAD",
       ifelse(grepl("^MSH0", Sample), "EEP", "SSP"))) %>%
  mutate(Sample_ID = case_when(Sample == "G449" ~ "G449",
                               Sample == "G445" ~ "G445",
                               Sample == "MSH009" ~ "#35552",
                               Sample == "MSH005" ~ "#34412",
                               Sample == "MSH638" ~ "#33556",
                               Sample == "MSH645" ~ "#36948"))

#~~ Plot for paper

colpal <- c(rev(wes_palette("IsleofDogs1", n = 5, type = "discrete")),
            wes_palette("GrandBudapest1", n = 4, type = "discrete"))
colpal <- colpal[c(1,2,3,5,4,6)]
cbPalette <- wes_palette("Moonrise2")

fac <- c("a","b","c")

psmc_df <- psmc_df %>%
  mutate(fac = as.factor(rep(fac, each = 762)))


fig_3 <- ggplot(mapping = aes(dplyr::filter(psmc_df, Sample == "MSH005" & mu == mid_mu & g == mid_g))) +
  geom_line(aes(YearsAgo, Ne, group = Boot, col = Sample), 
            dplyr::filter(boot_df, Sample == "MSH005" & mu == mid_mu & g == mid_g), 
            size = .05, alpha = 0.4, col = cbPalette[2]) +
  geom_line(aes(YearsAgo, Ne, col = fac), 
            dplyr::filter(psmc_df, Sample == "MSH005"), size = 0.7, alpha=0.6) +  
  theme_emily() +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 9)) +
  scale_colour_manual(values = cbPalette) +
  theme(legend.position = "none") +
  scale_x_log10(#breaks = trans_breaks("log10", function(x) 10^x),
               # labels = trans_format("log10", math_format(10^.x)),
                label = comma,
                breaks = c(10000,100000,1000000)) +
  scale_y_continuous(label = comma) +
  xlab("Years before present") + 
  ylab (expression(paste("IICR (scaled in units of 4",italic(N[e]),mu,")"))) +
  #ylab("Effective population size"~italic((N[e]))) +
  annotate("rect", xmin = 10000, xmax = 120000, ymin = -Inf, ymax = Inf,
           alpha = .1) +
  geom_vline(xintercept=22000, linetype = "dashed", colour = "grey") +
  annotate("text", x = 8700, y = 112000, label = "HOLOCONE", angle = 90, size = 2.5) +
  geom_vline(xintercept=2600000, colour = "grey") +
  annotate("text", x = 19000, y = 115000, label = "LGM 22 ka", angle = 90, size = 2.5) +
  geom_vline(xintercept=10000, colour = "grey") +
  annotate("text", x = 2250000, y = 25000, label = "PLEISTOCENE", angle = 90, size = 2.5) +
  geom_vline(xintercept=5300000, colour = "grey") +
  annotate("text", x = 4550000, y = 18500, label = "PLIOCENE", angle = 90, size = 2.5)


#ggsave("figs/Figure_3.tiff", fig_3, height = 4, width = 7)
ggsave("figs/Figure_3.tiff", fig_3, height = 6.5, width = 11, units = "cm")


# With legend

png(file="figs/psmc_figure_legend.png", units = "in", res = 300, height=5, width=9)

ggplot(mapping = aes(dplyr::filter(psmc_df, Sample == "MSH005" & mu == mid_mu & g == mid_g))) +
  geom_line(aes(YearsAgo, Ne, group = Boot, col = Sample), 
            dplyr::filter(boot_df, Sample == "MSH005" & mu == mid_mu & g == mid_g), 
            size = .1, alpha = 0.4, col = cbPalette[2]) +
  geom_line(aes(YearsAgo, Ne, col = fac), 
            dplyr::filter(psmc_df, Sample == "MSH005"), size = 0.8, alpha=0.6) +  
  theme_emily() +
  scale_colour_manual(values = cbPalette,
                      label = c(expression(paste("Low: ",mu, "= 0.8e-8, yr/gen = 3")), 
                                expression(paste("Medium: ",mu, "= 1.1e-8, yr/gen = 5.9")), 
                                expression(paste("High: ",mu, "= 1.3e-8, yr/gen = 10"))),
                      breaks = c("a", "b", "c"),
                      name = "Scaling") +
  theme(legend.position = "right") +
  scale_x_log10(#breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x)),
    label = comma,
    breaks = c(10000,100000,1000000)) +
  scale_y_continuous(label = comma) +
  xlab("Years before present") + 
  ylab (expression(paste("Inverse coalescence rate (scaled in units of 4",italic(N[e]),mu,")"))) +
  #ylab("Effective population size"~italic((N[e]))) +
  annotate("rect", xmin = 10000, xmax = 120000, ymin = -Inf, ymax = Inf,
           alpha = .1) +
  geom_vline(xintercept=22000, linetype = "dashed", colour = "grey") +
  annotate("text", x = 20000, y = 124000, label = "LGM 22 ka", angle = 90, size = 3) +
  geom_vline(xintercept=10000, colour = "grey") +
  annotate("text", x = 9300, y = 123000, label = "HOLOCONE", angle = 90, size = 3) +
  geom_vline(xintercept=2600000, colour = "grey") +
  annotate("text", x = 2400000, y = 121000, label = "PLEISTOCENE", angle = 90, size = 3) +
  geom_vline(xintercept=5300000, colour = "grey") +
  annotate("text", x = 4900000, y = 125000, label = "PLIOCENE", angle = 90, size = 3)

dev.off()


# Plot for supplementary

psmc_df <- psmc_df %>%
  mutate(pop = as.factor(pop),
         Sample_ID = factor(Sample_ID, levels = c("#33556", "#36948",
                                                     "#34412", "#35552",
                                                     "G445", "G449")))

boot_df <- boot_df %>%
  mutate(Sample_ID = factor(Sample_ID, levels = c("#33556", "#36948",
                                                  "#34412", "#35552",
                                                  "G445", "G449")))

png(file="figs/psmc_supp.png", units = "in", res = 300, height=8, width=8)

ggplot(mapping = aes(dplyr::filter(psmc_df, mu == mid_mu & g == mid_g))) +
  geom_line(aes(YearsAgo, Ne, col = pop), dplyr::filter(psmc_df, mu == mid_mu & g == mid_g), 
            size = 0.8, alpha=1) +  
  facet_wrap(~Sample_ID, ncol = 2) +
  geom_line(aes(YearsAgo, Ne, group = Boot, col = pop), boot_df, size = .1, alpha = 0.4) +
  theme_emily() +
  #theme(strip.text.x = element_blank()) +
  scale_colour_manual(values = cbPalette[c(1,2,3)],
                      name = "Population") +
  #theme(legend.position = "none") +
  scale_x_log10(#breaks = trans_breaks("log10", function(x) 10^x),
    # labels = trans_format("log10", math_format(10^.x)),
    label = comma,
    # limits = c(-1,5000000),
    breaks = c(10000,100000,1000000)) +
  scale_y_continuous(label = comma,
                     breaks = c(50000,100000)) +
  xlab("Years before present") + 
  ylab(expression(paste("IICR (scaled in units of 4",italic(N[e]),mu,")")))
  #ylab("Effective population size"~italic((N[e])))

dev.off()
