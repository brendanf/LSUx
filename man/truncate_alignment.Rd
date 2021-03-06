% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsux.R
\name{truncate_alignment}
\alias{truncate_alignment}
\title{Truncate an annotated Stockholm alignment file with secondary structure}
\usage{
truncate_alignment(alnfile, outfile, start = 1L, stop = 1000000L)
}
\arguments{
\item{alnfile}{(\code{character} filename or \code{\link{connection}}) the
stockholm alignment to be truncated.}

\item{outfile}{(\code{character} filename or \code{\link{connection}}) path
or connection to output the truncated alignment to.}

\item{start}{(\code{integer} scalar) first base to include in the truncated
alignment.}

\item{stop}{(\code{integer} scalar) last base to include in the truncated
alignment.}
}
\value{
NULL (invisibly)
}
\description{
Simply truncating a Stockholm alignment file that includes secondary
structure annotations may result in invalid base pairing, when only one half
of a pair is removed. This function converts these half-pairs to "X",
indicating unpairable bases.
}
\examples{
aln <- system.file(file.path("extdata", "fungi_32S.cm"), package = "LSUx")
truncate_alignment(aln, tempfile("trunc", fileext = ".stk"), 1, 500)
}
