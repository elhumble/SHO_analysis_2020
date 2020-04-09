library(data.table)
library(dplyr)
library(tidyr)
library(vcfR)
library(ape)

#~~ Summary stats

# Read in text file with output from samtools depth (manual)

stats <- read.table("data/mito/mito_nreads_focal.txt", header = T)

summary(stats$nreads)
summary(stats$coverage)

# Read in alignment of 6 focal individuals

t <- read.dna("data/mito/6_SHO_alignment.fasta", format = "fasta")
t <- ape::as.alignment(t)
t <- ape::as.DNAbin.alignment(t)

length(seg.sites(t))

dist.matrix <- ape::dist.dna(t, model = "raw", as.matrix = TRUE)

dist.df <- as.data.frame(dist.matrix)

sim_mat <- as.data.frame(dist.matrix) %>%
  gather() %>%
  mutate(similarity = (1-(value))*100) %>%
  mutate(ID1 = rownames(.))

summary(sim_mat)

