% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsux.R
\name{merge_5_8S}
\alias{merge_5_8S}
\title{Repair poor 5.8S detection in ITSx results using Infernal results.}
\usage{
merge_5_8S(itsx_result, csearch_result)
}
\arguments{
\item{itsx_result}{(\code{data.frame}) "positions" results from
\code{\link[rITSx:itsx]{rITSx::itsx()}}}

\item{csearch_result}{(\code{data.frame}) results from
\code{\link[inferrnal:cmsearch]{inferrnal::cmsearch()}} run on the same set of sequences,
using a 5.8S RNA model such as
\href{https://rfam.xfam.org/family/RF00002/cm}{RF00002}}
}
\value{
a \code{\link[tibble:tibble]{tibble::tibble}} in the same format output by
\code{\link[rITSx:itsx]{rITSx::itsx()}}, with updated positions for the boundaries of
5.8S.
}
\description{
The 3' end of the 5.8S RNA is quite variable, and it is sometimes not
detected by the HMM used in ITSx.
Using a CMs, Infernal is able to more reliably delimit 5.8S.
This function uses the results of \code{\link[inferrnal:cmsearch]{inferrnal::cmsearch()}} to
fill in missing 5.8S annotations in ITSx results.  It also updates the end of
ITS1 and beginning of ITS2 to match.
}
\examples{
# load sample data from inferrnal
seqfile <- system.file(
    file.path("extdata", "sample.fasta"),
    package = "inferrnal"
)
# ITSx has trouble with some of the reads
seq <- Biostrings::readDNAStringSet(seqfile)[c(1,4,20,32,33,43,46,49)]
# the result from ITSx is included as a dataset to avoid a package
# dependency, but this is the code to generate it.
#itsx_result <- rITSx::itsx(
#    in_file = seq,
#    positions = TRUE,
#    complement = FALSE,
#    cpu = 1,
#    read_function = Biostrings::readDNAStringSet
#)
pos <- itsx_result$positions
pos[pos$region == "5_8S",]
# find 5.8S using cmsearch
cm_5_8S <- system.file(
    file.path("extdata", "RF00002.cm"),
    package = "inferrnal"
)
cm_result <- inferrnal::cmsearch(seq, cm = cm_5_8S, cpu = 1)
# combine the results
merge_pos <- merge_5_8S(pos, cm_result)
merge_pos[merge_pos$region == "5_8S",]
}
