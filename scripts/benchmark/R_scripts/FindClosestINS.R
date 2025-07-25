#!/usr/bin/env Rscript

#
#
#
#
#
#

# (c) 2018-2021 SVUnit Authors. All Rights Reserved.
#
#
#
#
#
#
#
#
#
#

# This script is to identify the best matched SV from a reference callset (e.g. gnomAD) with a query SV callset based on bedtools closest results.
# Two SVs are considered matched if they have same SVTYPE and start coordinate difference < 500bp
# The script takes bedtools closest results (bed_a: query callset, bed_b: reference callset)
# The script will output a file with query SV id and the best matched reference SV id, and AF of the matched reference SV.

# a small value
MIN_VAL <- 1e-9

# get arguments
args <- commandArgs(trailingOnly=TRUE)

# input bed file
input_bed <- args[1]
# output file
output_file <- args[2]

# read bed file
dat <- read.table(input_bed, header=T, sep='\t', comment.char = '')
dat$dist <- abs(dat$start_a - dat$start_b)
dat$type_match <- ifelse(dat$svtype_a==dat$svtype_b, 1, 0)
dat$match <- ifelse(dat$dist<500 & dat$type_match==1, 1, 0)
dat <- dat[dat$match==1,]

# If a query SV is matched with multiple reference SVs, choose the one with smallest start distance.
# If there is a tie, choose the one with highest AF.
# If there is still a tie, choose the first one on the list.
dat <- dat[order(dat$query_svid, dat$dist, -dat$AF, decreasing = F),]
dat <- dat[!duplicated(dat$query_svid),]

dat <- dat[,c('query_svid','svid_b')]
write.table(dat, file=output_file, sep='\t', row.names=F, col.names=T, quote=F)


