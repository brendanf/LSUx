% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amplicon_extract.R
\name{find_primer}
\alias{find_primer}
\title{Find the best location for a primer sequence in a gappy alignment}
\usage{
find_primer(ungapped, gaps, primer)
}
\arguments{
\item{ungapped}{(\code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object)
The unaligned, gap-free sequences in the alignment.}

\item{gaps}{(\code{list} of \code{\link[IRanges:IRanges-class]{IRanges}} objects)
The locations of gaps in each sequence of \code{ungapped} when aligned.}

\item{primer}{(\code{\link[Biostrings:XString-class]{DNAString}} object)
The primer sequence to search for. It should be in the orientation which
matches the sequences in the alignment. (I.e., reverse primer sequences
should be reverse complemented prior to calling \code{find_primer()}.)}
}
\value{
named \code{integer} of length two, with elements \code{start} and \code{end} giving
the best fit range for the primer in the alignment.
}
\description{
This is a relatively simple algorithm which aligns the primer sequences to
each sequence in the alignment, maps the aligned locations into the alignment,
and takes the most frequent location. A warning is issued if more than 10\% of
the sequences in the alignment have non-consensus primer positions; a
frequent cause of this is if some of the sequences are fragmentary and do not
actually include the primer site, but in any case it is recommended to
manually check results.
}
