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
# Two SVs are considered matched if they have same SVTYPE, reciprocal overlap >50% and size similarity >50%
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

# calculate reciprocal overlap and size similarity
#dat$self_overlap_a <- (dat$end_a - dat$start_a)/dat$len_a
#dat$self_overlap_a[dat$len_a==0] <- 1
#dat$self_overlap_b <- (dat$end_b - dat$start_b)/dat$len_b
#dat$self_overlap_b[dat$len_b==0] <- 1
dat$ro <- dat$overlap_len/pmax(dat$len_a, dat$len_b)
dat$sim <- pmin(dat$len_a, dat$len_b)/pmax(dat$len_a, dat$len_b)

# SV type match
dat$type_match <- ifelse(dat$svtype_a==dat$svtype_b, 1, 0)
# SVs with reciprocal overlap > 50% and size similarity > 50% are considered matched
dat$match <- ifelse(dat$ro>0.5 & dat$sim>0.5 & dat$type_match==1, 1, 0)

# filter bed file
dat <- dat[dat$match==1,]

# distance between start and end of two SVs
dat$dist <- abs(dat$start_a - dat$start_b) + abs(dat$end_a - dat$end_b)

# If a query SV is matched with multiple reference SVs, choose the one with smallest start/end distance.
# If there is a tie, choose the one with highest AF.
# If there is still a tie, choose the first one on the list.
dat <- dat[order(dat$query_svid, dat$dist, -dat$AF, decreasing = F),]
dat <- dat[!duplicated(dat$query_svid),]

# output file with query svid and matched gnomad svid
dat <- dat[,c('query_svid','svid_b')]
write.table(dat, file=output_file, sep='\t', row.names=F, col.names=T, quote=F)


