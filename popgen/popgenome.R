setwd("/home/ubelix/ips/rchoudhury/Data/biscut_ind/")
library(PopGenome)
library(reshape)

vcf <- readData("popgen", format = "VCF",include.unknown=T)