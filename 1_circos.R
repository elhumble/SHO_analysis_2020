# Code for circos plotting

library(data.table)
library(dplyr)
options(scipen = 999)
library(RColorBrewer)
library(car)
source("scripts/RCircos.Long.Genome.Link.Plot.R")
library(RCircos)
library(purrr)

#~~~~~~~~~~~~
# Read in coordinate files parsed from LAST output (split by species, order retained)

cow <- fread("data/BosTau_coords.txt",
             col.names = c("Contig", "Start", "AlignLength", "Strand", "Length"))

sho <- fread("data/SHO_coords.txt", 
             col.names = c("Contig", "Start", "AlignLength", "Strand", "Length")) %>%
  mutate(Contig = gsub("_quiver_pilon", "", Contig))


#~~~~~~~~~~~~
# Get chromsome length stats

sho_align <- read.table("data/oryx_lengths", header = F) %>%  
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(sho, by = "Contig") %>%
  mutate(End = Start + AlignLength) 
#.[c(1,3,7)]


cow_align <- read.table("data/BosTau_lengths", header = F) %>%
  `colnames<-`(c("Contig", "Length")) %>%
  mutate(Contig = gsub(">", "", Contig)) %>%
  mutate(Length = as.numeric(Length)) %>%
  right_join(cow, by = "Contig") %>%
  mutate(End = Start + AlignLength)
#.[c(1,3,7)]


#~~~~~~~~~~~~
# Join both dataframes
# Filter by alignment length

align <- cbind(sho_align, cow_align) %>%
  .[c(1,2,3,4,7,8,9,10,14)] %>%
  `colnames<-`(c("Chromosome", "SHOContigLength", "chromStart", "AlignLength", "chromEnd", 
                 "Chromosome.1", "CowContigLength", "chromStart.1", "chromEnd.1"))

# get alignment summary stats

sum <- dplyr::select(align, Chromosome, SHOContigLength, AlignLength) %>%
  distinct(Chromosome, SHOContigLength, .keep_all = T) %>%
  dplyr::arrange(desc(SHOContigLength)) %>%
  top_n(29, SHOContigLength)

summary(sum$AlignLength)/1000
summary(sum$AlignLength)

# get SHO autosomes

top_chr <- dplyr::select(align, Chromosome, SHOContigLength) %>%
  distinct(Chromosome, SHOContigLength) %>%
  dplyr::arrange(desc(SHOContigLength)) %>%
  top_n(29, SHOContigLength)

summary(top_chr$SHOContigLength)

align <- top_chr %>%
  left_join(align, by = "Chromosome") %>%
  .[-c(2,3,8)] %>%
  mutate(Chromosome = gsub("Contig", "chr", Chromosome)) %>%
  filter(grepl("^GK", Chromosome.1)) %>%
  filter(AlignLength >= 10000) %>% # filter for alignment length
 # filter(AlignLength >= 5000) %>% # filter for alignment length
  dplyr::select(-AlignLength)


# Error checking
link.lengths <- align$chromEnd.1 - align$chromStart.1
which(link.lengths==0)


# Recode cow chromosome names
align$Chromosome.1 <- recode(align$Chromosome.1, '"GK000001.2" = "chr1"; "GK000002.2" = "chr2"; 
                             "GK000003.2" = "chr3"; "GK000004.2" = "chr4"; 
                             "GK000005.2" = "chr5"; "GK000006.2" = "chr6"; 
                             "GK000007.2" = "chr7"; "GK000008.2" = "chr8"; 
                             "GK000009.2" = "chr9"; "GK000010.2" = "chr10"; 
                             "GK000011.2" = "chr11"; "GK000012.2" = "chr12"; 
                             "GK000013.2" = "chr13"; "GK000014.2" = "chr14"; 
                             "GK000015.2" = "chr15"; "GK000016.2" = "chr16"; 
                             "GK000017.2" = "chr17"; "GK000018.2" = "chr18"; 
                             "GK000019.2" = "chr19"; "GK000020.2" = "chr20"; 
                             "GK000021.2" = "chr21"; "GK000022.2" = "chr22"; 
                             "GK000023.2" = "chr23"; "GK000024.2" = "chr24"; 
                             "GK000025.2" = "chr25"; "GK000026.2" = "chr26";
                             "GK000027.2" = "chr27"; "GK000028.2" = "chr28"; 
                             "GK000029.2" = "chr29"; "GK000030.2" = "chrX"')

# Recode SHO chromosome names
align <- align %>%
  mutate(Chromosome = gsub("HiC_scaffold_", "chr", Chromosome))


# Get link colours
color <- colorRampPalette(brewer.pal(10,"Paired"))(40)

# Wes
# color <- c("#9986A5", "#C6CDF7", "#7294D4", "#046C9A",  "#46ACC8","#5BBCD6", 
#            "#ABDDDE", "#00A08A", "#0B775E", "#81A88D", "#A2A475",
#            "#9A8822", "#CEAB07", "#E2D200", "#FAD510", "#F2AD00", 
#            "#D69C4E", "#F1BB7B", "#F5CDB4", "#D8A499", "#F8AFA8",  
#            "#F4B5BD", "#E6A0C4", "#FD6467", "#F2300F","#C93312","#9B110E", 
#            "#B24D17","#D67236", "#F98400")

# color <- sample(color) # get a good sample
pie(rep(1, length(color)), col = color , main="") 

# Or read in good color swatch
# color <- read.table("colors_1")
# color <- my_color$x

# Generate link.data
linkcolours <- matrix(nrow = 40,
                      ncol = 2)
linkcolours[,1] <- color
linkcolours[,2] <- c(paste0("chr", seq(1:38)), "chrX", "chrMT")
colnames(linkcolours) <- c("PlotColor", "Chromosome.1")
linkcolours <- as.data.frame(linkcolours)

link.data <- left_join(align, linkcolours, by = "Chromosome.1")

# Wes

# linkcolours <- matrix(nrow = 27,
#                       ncol = 2)
# linkcolours[,1] <- color
# linkcolours[,2] <- c(paste0("chr", seq(1:25)), "chrX", "chrMT")
# colnames(linkcolours) <- c("PlotColor", "Chromosome.1")
# linkcolours <- as.data.frame(linkcolours)
# 
# link.data <- left_join(align, linkcolours, by = "Chromosome.1")


#~~~~~~~~~~~~
# Make Ideograms

# SHO
sho_ideogram <- read.table("data/oryx_lengths", header = F) %>%
  `colnames<-`(c("Chromosome", "ChromEnd")) %>%
  mutate(ChromStart = 0) %>%
  dplyr::arrange(desc(ChromEnd)) %>%
  top_n(29, ChromEnd) %>%       
  .[c(1,3,2)] %>%
  mutate(Band = "qA1", Stain = c(rep(c("gneg", "gpos"),14), "gneg")) %>%
  mutate(Chromosome = gsub("HiC_scaffold_", "chr", Chromosome))

# reorder SHO chromosomes

sho_ideogram <- sho_ideogram[c(1,2,5,4,6,7,9,8,13,11,10,14,15,12,16,17,18,21,23,19,20,24,25,22,26,28,29,27,3),]


# Cow
cow_ideogram <- read.table("data/BosTau_lengths", header = F) %>%
  `colnames<-`(c("Chromosome", "ChromEnd")) %>%
  mutate(Chromosome = gsub(">", "", Chromosome)) %>%
  mutate(ChromStart = 0) %>%
  .[c(1,3,2)] %>%
  filter(grepl("^GK", Chromosome)) %>%
  mutate(Band = "qA1", Stain = rep(c("gneg", "gpos"),15))

cow_chr <- c(paste0("chr", seq(1:29)), "chrX")
cow_ideogram$Chromosome <- cow_chr

cow_ideogram <- map_df(cow_ideogram, rev)

# Correct species names
species.list <- c("SHO", "C")
cyto.list <- list(sho_ideogram, 
                  cow_ideogram)
RCircos.Multiple.Species.Core.Components(cyto.list, 
                                         species.list, NULL, 1, 0)

link.data[,1] <- paste(species.list[1], link.data[,1], sep="")
link.data[,4] <- paste(species.list[2], link.data[,4], sep="")




#~~~~~~~~~~~~
# Make Circos Plot

# Set plotting parameters
params <- RCircos.Get.Plot.Parameters()
params$base.per.unit <- 300000
params$highlight.width <- 0.7
params$line.color <- "gray"
RCircos.Reset.Plot.Parameters(params)

# Initialize graphic device (GUI or image file)
png(file="figs/circos.png", units = "in", res = 300, height=8, width=8)

par(cex=0.75)
RCircos.Set.Plot.Area()
# Plot chromosome ideogram and link lines
RCircos.Chromosome.Ideogram.Plot()
track.num <- 1
#RCircos.Link.Plot(link.data, track.num, FALSE)
RCircos.Long.Genome.Link.Plot(link.data, track.num = 1)

dev.off()

# 
# jpeg(file="figs/circos.jpg", units = "in", res = 300, height=8, width=8)
# 
# par(cex=0.75)
# RCircos.Set.Plot.Area()
# # Plot chromosome ideogram and link lines
# RCircos.Chromosome.Ideogram.Plot()
# track.num <- 1
# RCircos.Link.Plot(link.data, track.num, FALSE)
# #RCircos.Long.Genome.Link.Plot(link.data, track.num = 1)
# 
# dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#       Alignment Stats         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

detach(package:plyr)

# Using align df which includes NC dog chrosomes and longest 40 seal scaffs
# Contigs mapping to just one unique Dog chromsome

nchr <- align %>%
  dplyr::group_by(Chromosome) %>%
  dplyr::summarise(length(unique(Chromosome.1)))

nrow(filter(nchr, `length(unique(Chromosome.1))` == 1)) # 26
nrow(filter(nchr, `length(unique(Chromosome.1))` == 1)) / nrow(nchr)  # 65 %

# Contigs mapping mostly to one Dog chromosome

# number of hits to each dog chr
all_hits <- align %>%
  dplyr::group_by(Chromosome, Chromosome.1) %>%
  dplyr::summarise(n = length(Chromosome.1))

# number of total hits

total_hits <- all_hits %>%
  group_by(Chromosome) %>%
  dplyr::summarise(sum = sum(n))

# proportion of total hits for each dog chr

# largest number of hits

biggest_hit <- all_hits %>% 
  group_by(Chromosome) %>%
  filter(n == max(n))

df <- left_join(total_hits, biggest_hit, by = "Chromosome") %>%
  mutate(prop = n/sum)

nrow(filter(df, prop > 0.99))
nrow(filter(df, prop > 0.90))
nrow(filter(df, prop > 0.80))
nrow(filter(df, prop > 0.70))
nrow(filter(df, prop > 0.30))

