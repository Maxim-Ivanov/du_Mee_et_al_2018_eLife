# The function getOverlappingScores() integrates the sequencing
# coverage signal over each genomic interval of interest by
# calling the findOverlaps() function from GenomicRanges; 
# Input 1: intervals of interest (as GRanges object);
# Input 2: signal (as GrangesList; each GRanges corresponds to a
# separate track file, e.g. Bedgraph or Bigwig);
# The trim.names argument allows for stripping characters from
# the beginning and the end of each filename;
# Output: Granges object with additional metadata columns
# containing the integrated signal;

getOverlappingScores <- function(intervals, signal_grl, \
trim.names=c(0, 0)) {
  require(GenomicRanges)
  before <- trim.names[[1]]; after <- trim.names[[2]]
  for (i in seq_along(signal_grl)) {
    sample <- names(signal_grl)[[i]]
    sample <- substr(sample, before+1, nchar(sample)-after)
    cat(sample, "\n"); flush.console()
    signal <- signal_grl[[i]]
    hits <- findOverlaps(intervals, signal)
    query <- intervals[queryHits(hits)]
    subject <- signal[subjectHits(hits)]
    trimmed <- restrict(subject, start(query), end(query), \
keep.all.ranges=TRUE)
    by_obj <- by(trimmed, list(queryHits(hits)), \
function(x) { sum(width(x) * score(x)) })
    res <- t(do.call(rbind, list(by_obj)))
    out <- matrix(nrow=length(intervals), ncol=1, data=0)
    out[as.numeric(rownames(res))] <- res[, 1]
    mcols(intervals)[,sample] <- as.numeric(out)
  }
return(intervals)
}
