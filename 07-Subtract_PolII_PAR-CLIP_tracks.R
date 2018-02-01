# This script generates the difference track between
# Nrd1-AnchorAway and wild-type PolII PAR-CLIP samples
# (Schaughency et al., 2014 - accession GSE56435);

# Load the required packages:
library(rtracklayer)
library(GenomicRanges)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

# Load the PolII PAR-CLIP Bigwig files provided by authors
# (the original filenames were changes for convenience):

sq <- seqinfo(Scerevisiae)
aa_fw <- import.bw("Nrd1_AA_fw.bw", seqinfo = sq)
aa_rev <- import.bw("Nrd1_AA_rev.bw", seqinfo = sq)
neg_fw <- import.bw("Rpb2_withRap_fw.bw", seqinfo = sq)
neg_rev <- import.bw("Rpb2_withRap_rev.bw", seqinfo = sq)

# Subtract the negative control track from the Nrd1-AA track:
writeDeltaTrack <- function(gr1, gr2, fname, round.precision=3) {
  aa_rle <- coverage(gr1, weight="score")
  neg_rle <- coverage(gr2, weight="score")
  delta_rle <- round(aa_rle - neg_rle, round.precision)
  delta <- bindAsGRanges(score = delta_rle)
  con <- gzfile(fname, "w")
  writeLines("track type=bedGraph color=0,0,0, altColor=128,128,128", con)
  export.bedGraph(delta, con)
  close(con)
}
writeDeltaTrack(aa_fw, neg_fw, "Delta_AA_vs_negCtrl_fw.bg.gz")
writeDeltaTrack(aa_rev, neg_rev, "Delta_AA_vs_negCtrl_rev.bg.gz")

# Observe that the difference tracks are strand-specific.
# Positive and negative values in these tracks show the direction
# of difference between the Nrd1-AnchorAway and the control
# samples, not the strandness of PolII PAR-CLIP reads. Positive
# values mean that the level of nascent transcription was
# increased upon the nuclear depletion of Nrd1;