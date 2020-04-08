library(data.table)
library(dplyr)
library(tidyr)
library(vcfR)

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
  mutate(ID1 = rownames(dist.df))

summary(sim_mat)


#~~ VCF

vcf_file <- "data/mito/SHO_mito_hass_Q.recode.vcf.gz"
vcf <- read.vcfR(vcf_file, verbose = FALSE )

fix <- vcf@fix %>%
  as.data.frame() %>%
  dplyr::select(CHROM, POS, REF, ALT)

gt <- extract.gt(vcf) %>%
  as.data.frame() %>%
  mutate(CHROM = rownames(.)) %>%
  gather(Sample, Genotype, -CHROM) %>%
  separate(CHROM, c("CHROM", "POS"), sep = "_")
  
ad <- extract.gt(vcf, element = 'AD') %>%
  as.data.frame() %>%
  mutate(CHROM = rownames(.)) %>%
  gather(Sample, Depth, -CHROM) %>%
  separate(CHROM, c("CHROM", "POS"), sep = "_") %>%
 # left_join(fix, by = c("CHROM", "POS"))
  separate(Depth, c("Ref", "Alt")) %>%
  mutate(Ref = as.numeric(Ref),
         Alt = as.numeric(Alt)) %>%
  mutate(Call = case_when(Ref >= Alt ~ 0, 
                          TRUE ~ 1)) %>%
  left_join(gt, by = c("CHROM", "POS", "Sample")) %>%
  mutate(mismatch = case_when(Call != Genotype ~ "FLAG",
                              TRUE ~ "OK"))

depth <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)

# Haplotype network

library(ape)
library(pegas)

t <- read.dna("data/mito/All_CR_alignment.fasta", 
              format = "fasta")

ids <- labels(t)
names(t) <- ids

#source <- c(rep("SN", 3), rep("IY", 40),
#            "EAD","EAD","EEP",
 #           "EEP","SSP","SSP")

#labels(t) <- gsub("Y_.+", "Y", labels(t))

hap <- haplotype(t)
hap
hapnet <- haploNet(hap)

lab <- as.character(attr(hapnet, "freq"))
attr(hapnet, "labels") <- lab

ind.hap<-with(
  utils::stack(setNames(attr(hap, "index"), rownames(hap))),
  table(hap=ind, pop=names(t)[values]))

ind.hap <- as.data.frame(ind.hap) %>%
  filter(hap == "XXXIV")

plot(hapnet, size=attr(hapnet, "freq"), scale.ratio=0.4,
     labels = F, show.mutation = 0)

# StrataG

library(strataG)

cr <- read.fasta("data/mito/six_CR_sequences.fasta")

nucleotideDiversity(cr)

# tree

tree <- read.tree("data/mito/RAxML_bestTree.test_tree")
plot(tree)

data(dolph.seqs)
data(dolph.strata)
dloop <- df2gtypes(dolph.strata[, c("id", "fine", "id")], ploidy = 1,
                   schemes = dolph.strata[, c("fine", "broad")], sequences = dolph.seqs)
dloop <- labelHaplotypes(dloop, "Hap.")$gtypes

nucleotideDivergence(dloop)

