# This R pipeline was used to choose the candidate CUTs (n=68);

# Load the required package:
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Load all CUTs (n=925) from Xu et al., 2009 (Supplementary
# Table 3). Observe that the genomic coordinates were converted
# from sacCer2 to sacCer3 using kentUtils liftOver:
cuts <- read.table("Xu2009_CUTs_sacCer3.bed", sep="\t", header=F, stringsAsFactors=F)
cuts <- GRanges(seqnames=cuts$V1, ranges=IRanges(cuts$V2, end=cuts$V3, names=cuts$V4), strand=cuts$V6, seqinfo=seqinfo(Scerevisiae))

# Load all SUTs (n=847) from the same source:
suts <- read.table("Xu2009_SUTs_sacCer3.bed", sep="\t", header=F, stringsAsFactors=F)
suts <- GRanges(seqnames=suts$V1, ranges=IRanges(suts$V2, end=suts$V3, names=suts$V4), strand=suts$V6, seqinfo=seqinfo(Scerevisiae))

# Load all SGD genes (n=6692) from UCSC:
sgd <- read.table("sgdGene_UCSC.bed", sep="\t", header=F, comment.char="#", stringsAsFactors=F)
sgd <- GRanges(seqnames=sgd$V1, ranges=IRanges(sgd$V2, end=sgd$V3, names=sgd$V4), strand=sgd$V6, seqinfo=seqinfo(Scerevisiae))

# Combine CUTs, SUTs and SGD genes together:
all_genes <- c(cuts, suts, sgd)

# Remove CUTs which overlap any SUT or SGD gene on either strand:
cuts <- cuts[countOverlaps(cuts, all_genes, ignore.strand=T)==1]

# For each CUT, find the nearest downstream SGD gene on the same
# strand:
downstr <- sgd[precede(cuts, sgd)]

# Get the intervals between each CUTs and its downstream gene:
gaps <- pgap(cuts, downstr)

# Find CUT/gene pairs which meet the following conditions:
# i) (100 bp) <= (gap width) <= (1500 bp);
# No other transcription unit (CUT, SUT or SGD gene) overlapping
# the gap on either strand;
valid <- width(gaps)<=1500 & width(gaps)>=100 & countOverlaps(gaps, all_genes, ignore.strand=T)==0

# Subset CUTs, gaps and downstream genes by the logical index:
cuts <- cuts[valid]
gaps <- gaps[valid]
downstr <- downstr[valid]

# Save the GRanges objects for future use:
saveRDS(cuts, "cuts.RData")
saveRDS(gaps, "gaps.RData")
saveRDS(downstr, "downstr.RData")
