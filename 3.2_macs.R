library(ggplot2)
source("scripts/plot_psmc.R")
library(data.table)
library(plyr)
library(tidyr)
library(wesanderson)
library(scales)
options(scipen=999)
library(glue)
library(purrr)

# Write scripts for running macs simulations

#~~ Specify variables for plotting

i.iteration=20 # Number of iterations to use in file
s=100 # bin size

#~~ Read in main PSMC files

psmc_files <- paste("data/psmc/", list.files(path = "data/psmc", pattern="*.psmc"), sep = "")
dataset_names <- list.files(path = "data/psmc/", pattern="*.psmc")
dataset_names <-  gsub(".psmc", "", dataset_names)
names(psmc_files) <- dataset_names

#~~ Run for multiple values of mu and g:

mid_mu <- 1.1e-8
mid_g <- 6.2

psmc <- map_df(psmc_files, psmc.result, i.iteration, mu = mid_mu, s, g = mid_g, .id = "Sample") %>%
    mutate(pop = ifelse(grepl("^G", Sample) | grepl("^MSH2", Sample), "EAD",
                      ifelse(grepl("^MSH0", Sample), "EEP", "SSP")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#          Simulations           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Code modified from Annabel Beichman otter MBE paper

#~~ Function to write bash script

macs_script <- function(df, trim, nsims, len) {
  
  Sample <- unique(df$Sample)
  
  Ne <- df %>%
    dplyr::select(Ne)
  
  # Ne values until ancient time
  # scale it relative to oldest time in inference (could also do most recent size, just be consistent)
  # this is now relative to ancient size, so is ready for macs

  # Generation times: scale

  # if first time is -0, make it just 0.
  # times_gen_trimancient_4Na[0]=0 # set first time to zero if its -0.
  
  Na <- as.numeric(df[trim,3])
  
  scale <- df %>%
    slice(1:trim) %>%
    mutate(Ne_trimancient_Na = Ne/Na,
           times_gen = YearsAgo / g,
           times_gen_trimancient_4Na = times_gen / (4 * Na))
  
  r <- 1e-08 # recombination rate
  
  ss = 2 # two haplotypes (one genome)
  theta = 4 * Na * mid_mu 
  rho = 4 * Na * r 
  
  # Command
  
  rt <- "1:00:00"
  vmem <- "1G"
  array <- glue("1-{nsims}")
  outdir <-
    "/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2020/data/psmc/macs/"
  macs_path <-
    "/exports/cmvm/eddie/eb/groups/ogden_grp/software/macs-master/"
  sge <- "${SGE_TASK_ID}"
  
  # Use glue to stick everything together
  
  geo <-
    glue(
      "#!/bin/sh \n# Grid Engine Options \n#$ -N Macs \n#$ -cwd \n#$ -l h_rt={rt}",
      "\n#$ -l h_vmem={vmem} \n#$ -t {array} \n#$ -R y \n\n# Jobscript to sim chr with macs",
      "\n# Sample {Sample}",
      "\n. /etc/profile.d/modules.sh"
    )
  
  pre <-
    glue("{macs_path}macs {ss} {len} -t {theta} -r {rho} -s $SGE_TASK_ID") # -s $SGE_TASK_ID
  
  eN <-
    glue("-eN {scale$times_gen_trimancient_4Na} {scale$Ne_trimancient_Na}", sep = "")
  eN <- paste(eN, collapse = " ")
  out <- glue("> {outdir}macs_{trim}_sim_{sge}_{Sample}.txt")
  
  write.table(
    glue("{geo} \n\n{pre} {eN} {out}"),
    glue("macs_scripts/4.3_macs_{trim}_{Sample}.sh"),
    quote = F,
    col.names = F,
    row.names = F
  )
  
}

#purrr::walk(1:28, macs_script, sim[[1]], 100, 2500000) # walk for out file otherwise map

system("mkdir macs_scripts")

sim <- psmc %>%
  #dplyr::filter(Sample == "MSH005") %>%
  dplyr::group_by(Sample) %>%
  dplyr::distinct(Ne, YearsAgo, .keep_all= TRUE) %>% # remove duplicate Nes and times: check in plotting as well
  dplyr::distinct(Ne, .keep_all = T) 

sim %>%
  split(.$Sample) %>%
  map(macs_script, 28, 1000, 2500000)

# Transfer to server in terminal

# scp macs_scripts/4.3_macs* ehumble@eddie3.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2020/

# Run on Eddie

# for i in 4.3_macs_28_*; do qsub $i; done

# Copy out to local

# scp ehumble@eddie.ecdf.ed.ac.uk:/exports/cmvm/eddie/eb/groups/ogden_grp/emily/SHO_reseq_2020/data/psmc/macs/macs_* data/macs/


