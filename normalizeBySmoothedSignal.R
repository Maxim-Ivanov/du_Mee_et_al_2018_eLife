# This function is used to normalize one track (gr1) by smoothed
# values from the other track (gr2). To avoid division by zero,
# a small pseudocount is added to gr2. The smoothing is done by
# running mean with adjustable window size;
# The output data is also normalized to 1 million tags by
# default;
# The shift.start argument is used to avoid zero-width intervals;

normalizeBySmoothedSignal <- function(gr1, gr2, pscount=0.1, window=51, shift.start=1, round.precision=3, tpm.norm=TRUE) {
  rle1 <- mcolAsRleList(gr1, "score")
  rle2 <- mcolAsRleList(gr2, "score")
  rle2[is.na(rle2)] <- 0
  rle2 <- rle2 + pscount
  smooth_rle2 <- runmean(rle2, k=window, endrule="constant")
  scaled <- rle1 / smooth_rle2
  if (isTRUE(tpm.norm)) {
    total <- sum(unlist(scaled), na.rm=T) / 1000000
    scaled <- scaled / total
  }
  scaled <- round(scaled, round.precision)
  gr_out <- bindAsGRanges(score=scaled)
  strand(gr_out) <- strand(gr1)
  start(gr_out) <- start(gr_out) - shift.start
  seqinfo(gr_out) <- seqinfo(gr1)
  return(gr_out)
}
