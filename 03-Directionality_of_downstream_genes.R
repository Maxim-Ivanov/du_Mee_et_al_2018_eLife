# This R pipeline was used to assess the directionality of
# nascent transcription at genes which are located downstream
# from the candidate CUTs (n=68);

# Load the required packages:
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Load the required custom functions:
source("batchReadTrackData.R")
source("getOverlappingScores.R")

# Load the GRanges objects (were produced by the "01-Candidate_CUTs.R" pipeline):
cuts <- readRDS("cuts.RData")
gaps <- readRDS("gaps.RData")
downstr <- readRDS("downstr.RData")

# Shrink gaps and downstream genes to 100 bp (to consider nascent
# transcription only within 100 bp from the TSS):
dnc <- resize(gaps, 100, fix="end")
coding <- resize(downstr, 100)

# Flip the DNC (divergent non-coding transcription) regions to
# the opposite strand:
strand(dnc) <- ifelse(strand(dnc)=="+", "-", "+")

# Load names of Bedgraph files (were produced by the
# "02-Remapping_NET-Seq_raw_data.sh" pipeline):
bg_filenames <- list.files(".", pattern="norm1M.bg.gz$")

# Read the content of Bedgraph files to a GRangesList object:
bg_files <- batchReadTrackData(bg_filenames, \
seqinfo=seqinfo(Scerevisiae))

# Calculate the NET-Seq coverage over the 100 bp intervals around
# TSS of the downstream genes:
dnc_cov <- getOverlappingScores(dnc, bg_files)
coding_cov <- getOverlappingScores(coding, bg_files)

# Calculate the directionality index for the downstream genes;
# Observe that we sum the NET-Seq coverage values over all
# columns because the input Bedgraph files are replicates of the
# same sample;
# Also we add a small pseudocount to the DNC region to avoid
# division by zero:
dir_index <- rowSums(mcols(coding_cov)) / \
(rowSums(mcols(dnc_cov)) + 0.25)

# Save all results in a Tab-delimited file:
cuts_df <- as.data.frame(cuts)
coding_df <- as.data.frame(coding)
coding_df$gene_name <- rownames(coding_df)
dnc_df <- as.data.frame(dnc)
dist <- as.data.frame(width(gaps))
dir_index <- as.data.frame(dir_index)
output_df <- cbind(cuts_df, dist, dir_index, coding_df, dnc_df)
write.table(output_df, \
"Directionality_of_CUT_downstream_genes.txt", sep="\t", quote=F)

