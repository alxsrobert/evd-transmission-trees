if(!is.element("main_script.R", dir())) stop("Set working directory to file location") 

#### Import libraries and functions ####

library(igraph)
library(ape)
library(adegenet)
library(data.table)
library(dplyr)
library(outbreaker2)
library(ggplot2)

source("function_figures.R")
source("function_generate_outbreaks.R")

#### Define parameters ####
# Create a reference sequence of a,t,c, and g base
seq_ref <- sample(x = c("a", "t", "c", "g"), size = 18000, replace = T)
# Proportion of the sequence that can vary
prop_var <- .002
# Maximum number of imports
n_import <- 50
# time span of importations
date_range <- as.Date(c("2021-01-01", "2022-01-01"))
# Maximum number of cases
max_cases <- 300
# Parameters of the serial interval
mean_w <- 14.2
sd_w <- 7.1
shape_w <- (mean_w**2) / (sd_w**2)
scale_w <- (sd_w**2) / mean_w
w_dens <- pgamma(0:100, shape = shape_w, scale = scale_w) %>% diff
# Parameters of the incubation period
mean_f <- 9.1
sd_f <- 7.3
shape_f <- (mean_f**2) / (sd_f**2)
scale_f <- (sd_f**2) / mean_f
f_dens <- pgamma(0:100, shape = shape_f, scale = scale_f) %>% diff
# Mutation rate
mut_rate <- 0.2
# Reproduction number
r0 <- 1
# Report rate of the genetic sequences
prob_seq_rep <- .4

#### Generate outbreak ####

set.seed(1)
# Simulate the outbreaks, outbreaks contains the epi dataset, and the sequences
outbreaks <- generate_outbreak(n_import, date_range, max_cases, seq_ref, prop_var,
                          w_dens, f_dens, mut_rate, r0, prob_seq_rep)
# Select the chains with more than 12 genetic sequences reported 
chains <- which(table(as.numeric(outbreaks[[1]]$chain)) > 12)
# Keep the cases from "chains"
case_keep <- which(is.element(outbreaks[[1]]$chain, chains))
outbreaks[[1]] <- outbreaks[[1]][case_keep,]
outbreaks[[2]] <- outbreaks[[2]][case_keep,]

#### Save simulated sequences as fasta file ####

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    if(data[rowNum, "seq"] != ""){
      fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"ID"], sep = "")))
      fastaLines = c(fastaLines, as.character(data[rowNum,"seq"]))
    }
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

writeFasta(outbreaks[[2]], "example.fasta")

#### Import to analyse with outbreaker2 ####

sequences_fas <- read.FASTA("example.fasta")
lab_seq <- labels(sequences_fas)

#### Analysis with outbreaker2 ####
# Extract the columns ID, onset date, and chain number from outbreaks[[1]] 
combined_data <- outbreaks[[1]][is.element(ID, lab_seq),.(ID, date_ons, chain)]
combined_data <- combined_data[order(chain, date_ons),]
# Match the cases' ID in combined_data and sequence labels
seq_in_dt <- match(combined_data$ID, lab_seq)
seq_o2 <- sequences_fas[seq_in_dt]
names(seq_o2) <- seq_along(seq_o2)
# Create outbreaker_data object
data <- outbreaker_data(dna=seq_o2,
                        dates=as.numeric(combined_data$date_ons-
                                           min(combined_data$date_ons)),
                        w_dens = w_dens,
                        f_dens = f_dens,
                        id_in_dna = seq_along(seq_o2))
# Create outbreaker_config object
config <- create_config(n_iter= 5e4, data = data
                        , init_mu = 2e-5
                        , outlier_threshold = 3
                        , init_tree = "star"
)
library(tictoc)
# Run combined model with less stringent importation threshold
Ebola_out_o2_3 <- outbreaker(data = data, config = config)
# Run combined model with more stringent importation threshold
config$outlier_threshold <- 2
config <- create_config(config)
Ebola_out_o2_2 <- outbreaker(data = data, config = config)

# Run epi model with more stringent importation threshold
data$dna <- NULL
data$id_in_dna <- NULL
data <- outbreaker_data(data = data)
config$outlier_threshold <- 1.1
config <- create_config(config, data = data)
Ebola_out_WT_1 <- outbreaker(data = data, config = config)
# Run epi model with less stringent importation threshold
config$outlier_threshold <- 3
config <- create_config(config)
Ebola_out_WT_3 <- outbreaker(data = data, config = config)

#### Generate figures ####

## Time series to be generated using a report rate
figure1B(outbreaks)

## Alpha plot:
figure3(out_wt = Ebola_out_WT_3, out_comb = Ebola_out_o2_3, data = outbreaks, burnin = 2e3)

## Consensus trees:
figure4(out_wt = Ebola_out_WT_3, out_comb = Ebola_out_o2_3, data = outbreaks, burnin = 2e3)

## Alpha plots by chain
figure5(out_wt = Ebola_out_WT_3, out_comb = Ebola_out_o2_3, data = outbreaks, burnin = 2e3)

## Proportion linked to the "right" ancestor
figure6(out_wt = Ebola_out_WT_3, out_comb = Ebola_out_o2_3, data = outbreaks, burnin = 2e3)

## Supplementary Figure 2
supp_fig2(out_wt = Ebola_out_WT_3, out_comb = Ebola_out_o2_3, burnin = 2e3)

## Supplementary Figure 4
figure3(out_wt = Ebola_out_WT_1, out_comb = Ebola_out_o2_2, data = outbreaks, burnin = 2e3)

## Supplementary Figure 5
figure4(out_wt = Ebola_out_WT_1, out_comb = Ebola_out_o2_2, data = outbreaks, burnin = 2e3)

## Supplementary Figure 6
figure6(out_wt = Ebola_out_WT_1, out_comb = Ebola_out_o2_2, data = outbreaks, burnin = 2e3)
