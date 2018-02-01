# The Nrd1 PAR-CLIP signal (Schulz et al., 2013) is expected to
# depend on transcription rate of the underlying locus. The
# observed Nrd1 occupancy may be overestimated in highly
# transcribed genes and underestimated in low transcribed ones.
# To analyze the Nrd1 occulancy of CUT genes in an unbiased
# manner, the Nrd1 PAR-CLIP signal has to be normalized by the
# 4sU-Seq signal from the same samples.

# Load the required packages:
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Load the required custom functions:
source("batchReadTrackData.R")
source("getOverlappingScores.R")
source("")

# Load the unscaled Bedgraph file representing merged replicates
# of wild-type 4sU-Seq sample (this file was produced by the
# "05-Reanalysis_of_4sU-Seq_and_Nrd1_PAR-CLIP_raw_data.sh"
# pipeline):
thio_data <- import.bedGraph(4sU_wt_merged.bg.gz, seqinfo=seqinfo(Scerevisiae))
strand(thio_data) <- ifelse(score(thio_data) >= 0, "+", "-")
score(thio_data) <- abs(score(thio_data))

# Load the Nrd1 PAR-CLIP Bedgraph files which have to be
# normalized by the 4sU-Seq track:
parclip_filenames <- list.files(".", pattern="^PAR-CLIP_signal.*norm1M.bg.gz")
parclip_data <- batchReadTrackData(parclip_filenames, seqinfo=seqinfo(Scerevisiae))

# Save the original headers from PAR-CLIP Bedgraph files:
getHeader <- function(x) {
    con <- file(x, "r")
    header <- readLines(con, n=1)
    close(con)
    return(header)
}
parclip_headers <- lapply(parclip_filenames, getHeader)

# Now do the normalization for each Nrd1 PAR-CLIP sample. Observe
# that the mendoapply() function can be very slow:
parclip_data_norm <- mendoapply(normalizeBySmoothedSignal, parclip_data, thio_data)

# Save the resultant GRangesList object as gzipped Bedgraph
# files:
writeBgWithCustomHeader <- function(gr, fname, header) {
    df <- subset(as.data.frame(gr), select=c("seqnames", "start", "end", "score"))
    extension <- paste(".", tools::file_ext(name))
    fname <- gsub(extension, "_norm4sU_1M.bg.gz", name)
    con <- file(fname, "w")
    writeLines(header, con)
    write.table(df, con, sep="\t", quote=F, row.names=F, col.names=F)
    close(con)
}
mapply(writeBgWithCustomHeader, parclip_data_norm, names(parclip_data_norm), parclip_headers)

# Finally quantify the 4sU-/TPM-normalized Nrd1 PAR-CLIP
# occupancy of the candidate CUT genes (n = 68):
cuts <- readRDS("cuts.RData")
cuts_cov <- getOverlappingScores(cuts, parclip_data_norm)

# Save the results as Tab-delimited file:
write.table(as.data.frame(cuts_cov), file="Nrd1_PAR-CLIP_occupancy_of_68_CUTs.txt", sep="\t", quote=F)
