# This R pipeline was used to assess the Nrd1-dependent
# transcription termination capacity of the candidate CUTs
# (n=68).
# The PolII PAR-CLIP data from wild-type and Nrd1-AnchorAway
# samples (Schaughency et al., 2014 - accession GSE56435) were
# quantified over the whole CUT genes and their "gap" regions
# (the interval between CUT and the downstream gene on the same
# strand);

# Load the required packages:
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Load the required custom functions:
source("batchReadTrackData.R")
source("getOverlappingScores.R")

# Load the GRanges objects (were produced by the
# "01-Candidate_CUTs.R" pipeline):
cuts <- readRDS("cuts.RData")
gaps <- readRDS("gaps.RData")
downstr <- readRDS("downstr.RData")

# Load names of the original Bigwig files (downloaded from
# GSE56435). Observe that these files are strand-specific:
bw_filenames_fw <- list.files(".", pattern="fw.bw$")
bw_filenames_rev <- list.files(".", pattern="rev.bw$")

# Read the content of stranded Bigwig files to GRangesList
# objects:
bw_files_fw <- batchReadTrackData(bw_filenames_fw, strand="+", \
seqinfo=seqinfo(Scerevisiae))
bw_files_rev <- batchReadTrackData(bw_filenames_rev, \
strand="-", seqinfo=seqinfo(Scerevisiae))

# Combine stranded Bigwig files corresponding to the same sample:
bw_files <- GRangesList()
for (i in seq_along(bw_files_fw)) {
  fw <- bw_files_fw[[i]]
  rev <- bw_files_rev[[i]]
  name <- gsub("_fw.bw", "", names(bw_files_fw)[[i]])
  cat(name, "\n"); flush.console()
  bw_files[[name]] <- sort(c(fw, rev))
}

# Calculate the PolII PAR-CLIP coverage:
cuts_cov <- getOverlappingScores(cuts, bw_files)
gaps_cov <- getOverlappingScores(gaps, bw_files)

# For each PAR-CLIP sample (wild-type of Nrd1-AA) calculate the
# gap/CUT ratio (the readthrough index). Observe that custom
# column names were used (change if necessary):
norm_gaps_wt <- mcols(gaps_cov)$Rpb2_withRap / width(gaps_cov)
norm_gaps_AA <- mcols(gaps_cov)$Nrd1_AA / width(gaps_cov)
norm_cuts_wt <- mcols(cuts_cov)$Rpb2_withRap / width(cuts_cov)
norm_cuts_AA <- mcols(cuts_cov)$Nrd1_AA / width(cuts_cov)
readthrough_wt <- norm_gaps_wt / norm_cuts_wt
readthrough_AA <- norm_gaps_AA / norm_cuts_AA

# Finally calculate the increase of transcription readthrough
# upon Nrd1 depletion:
nrd1_dependency <- readthrough_AA / readthrough_wt

# Save the results in a Tab-delimited file:
output_df <- cbind(as.data.frame(cuts_cov), \
as.data.frame(gaps_cov), as.data.frame(nrd1_dependency))
write.table(output_df, "Nrd1_dependency_of_CUTs.txt", \
sep="\t", quote=F)
