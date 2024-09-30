#' Convert range from ungapped sequence to gapped sequence
#'
#' @param rstart
#' (`integer` scalar)
#' the start point of the range
#' @param rend
#' (`integer` scalar)
#' the end point of the range
#' @param gaps
#' ([`IRanges`][IRanges::IRanges-constructor] object)
#' the locations of gaps in the gapped sequence.
#'
#' @return named `integer` of length two, with elements `"start"` and `"end"`,
#' giving the coordinates in the gapped sequence which correspond to coordinates
#' `"rstart"` and `"rend"` in the ungapped sequence.
#' @export
gap_fill <- function(rstart, rend, gaps) {
  i = 1L
  wsum = 0L
  preInsert <- 0L
  postInsert <- 0L
  for (i in seq_along(gaps)) {
    wsum <- wsum + gaps@width[i]
    insert <- gaps@start[i] + gaps@width[i] - 1L - wsum
    if (insert < rstart) {
      preInsert <- wsum
    }
    if (insert < rend) {
      postInsert <- wsum
    } else {
      break
    }
  }
  c(
    start = preInsert + rstart,
    end = postInsert + rend
  )
}

#' Find the best location for a primer sequence in a gappy alignment
#' 
#' This is a relatively simple algorithm which aligns the primer sequences to
#' each sequence in the alignment, maps the aligned locations into the alignment,
#' and takes the most frequent location. A warning is issued if more than 10% of
#' the sequences in the alignment have non-consensus primer positions; a
#' frequent cause of this is if some of the sequences are fragmentary and do not
#' actually include the primer site, but in any case it is recommended to
#' manually check results.
#'
#' @param ungapped
#' ([`DNAStringSet`][Biostrings::XStringSet-class] object)
#' The unaligned, gap-free sequences in the alignment.
#' @param gaps
#' (`list` of [`IRanges`][IRanges::IRanges-class] objects)
#' The locations of gaps in each sequence of `ungapped` when aligned.
#' @param primer ([`DNAString`][Biostrings::XString-class] object)
#' The primer sequence to search for. It should be in the orientation which
#' matches the sequences in the alignment. (I.e., reverse primer sequences
#' should be reverse complemented prior to calling `find_primer()`.)
#'
#' @return named `integer` of length two, with elements `start` and `end` giving
#' the best fit range for the primer in the alignment.
#' @export
find_primer <- function(ungapped, gaps, primer) {
  aln <- Biostrings::pairwiseAlignment(
    ungapped,
    primer,
    type = "local-global",
    substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix()
  )
  result <- mapply(
    gap_fill,
    aln@pattern@range@start,
    aln@pattern@range@start + aln@pattern@range@width - 1L,
    gaps
  )
  start <- result[1,]
  end <- result[2,]
  bestscore <- max(aln@score)
  beststart <- as.integer(names(which.max(table(start))))
  bestend <- as.integer(names(which.max(table(end))))
  startmismatch <- start != beststart
  endmismatch <- end != bestend
  if (sum(startmismatch | endmismatch) > 0.1 * length(start))
    warning("more than 10% of sequences in reference alignment have variant\n",
            "locations for primer ", as.character(primer))
  c(start = beststart, end = bestend, score = bestscore)
}

as_StockholmMSA <- function(aln, name = "aln") {
  if (methods::is(aln, "StockholmDNAMultipleAlignment") ||
      methods::is(aln, "StockholmRNAMultipleAlignment")) {
    # ok, no problem
    aln
  } else if (methods::is(aln, "connection")) {
    assertthat::assert_that(summary(aln)[["can read"]] == "yes")
    inferrnal::read_stockholm_msa(aln)
  } else if (is.character(aln) &&
             length(aln) == 1 && 
             assertthat::is.readable(aln)) {
    inferrnal::read_stockholm_msa(aln)
  } else if (is.character(aln)) {
    # alignment given as character string
    assertthat::assert_that(
      dplyr::n_distinct(nchar(aln)) == 1
    )
    tryCatch(
      inferrnal::StockholmRNAMultipleAlignment(aln),
      error = function(e) inferrnal::StockholmRNAMultipleAlignment(aln)
    )
  } else if (methods::is(aln, "DNAStringSet") ||
             methods::is(aln, "DNAMultipleAlignment")) {
    inferrnal::StockholmDNAMultipleAlignment(aln)
  } else if (methods::is(aln, "RNAStringSet") ||
             methods::is(aln, "RNAMultipleAlignment")) {
    inferrnal::StockholmRNAMultipleAlignment(aln)
  } else {
    stop("'", name, "' should be a connection, a filename, or a DNA or RNA alignment")
  }
}

mark_ref_line <- function(aln, start, end, mark_char) {
    # ensure there is a valid reference line in the alignment
    if (!"RF" %in% names(aln@GC)) {
        # make a reference line out of the alignment consensus
        aln@GC$RF <- chartr(".", "-", aln@unmasked) |>
            Biostrings::consensusString() |>
            chartr(old = "-", new = ".")
    } else {
        if (grepl(mark_char, aln@GC$RF, fixed = TRUE)) {
            warning("warning:reference line of supplied alignment includes",
                    "character", shQuote(mark_char), ".")
        }
    }
    # find non-gap positions in the reference line
    refpos <- IRanges::gaps(
        Biostrings::matchPattern(".", aln@GC$RF),
        start = 1,
        end = Biostrings::nchar(aln@GC$RF)
    )
    
    # find the non-gap positions which match the primers
    mark_refpos <- IRanges::findOverlapPairs(
        refpos,
        IRanges::IRanges(start = start, end = end)
    )
    mark_refpos <- IRanges::pintersect(mark_refpos)
    
    aln@GC$RF <- Biostrings::replaceAt(
        aln@GC$RF,
        mark_refpos,
        c(strrep(mark_char, mark_refpos@width))
    )
    aln
}

#' Find the location of an amplicon in an alignment
#'
#' @param aln
#' (`connection`, `character` string giving a file name in Stockholm format,
#' [`DNAMultipleAlignment`][Biostrings::MultipleAlignment-class],
#' [`RNAMultipleAlignment`][Biostrings::MultipleAlignment-class],
#' [`DNAStringSet`][Biostrings::XStringSet-class],
#' [`RNAStringSet`][Biostrings::XStringSet-class],
#' [`StockholmMultipleAlignment`][inferrnal::StockholmMultipleAlignment-class],
#' or `character` vector)
#' DNA or RNA multiple alignment in which to search for an amplicon.
#' @param fwd_primer
#' (`character` string or [`DNAString`][Biostrings::XString-class])
#' Forward primer sequence to define the target amplicon.
#' @param rev_primer
#' (`character` string or [`DNAString`][Biostrings::XString-class])
#' Reverse primer sequence to define the target amplicon. Should be given from
#' 5' to 3' in the primer; i.e. the reverse complement of the expected sequence
#' in the alignment.
#' @param trim
#' (one of `"none"`, `"retain"`, or `"remove"`)
#' Choice of how to trim the alignment: if `"none"` then the alignment is not
#' trimmed; if `"retain"` then the alignment is trimmed to the amplicon,
#' including the primer sites; if `"remove"` then the alignment is trimmed to
#' the amplicon and the primer sites are also removed.
#' @param mark
#' (`logical` flag)
#' If `TRUE` (default) the primer sites are marked in the alignment RF line.
#' @param fwd_char
#' (single `character`)
#' Character to use for marking the forward primer location in the RF line.
#' @param rev_char
#' (single `character`)
#' Character to use for marking the reverse primer location in the RF line.
#' @param outfile
#' (`character` file name or [`connection`])
#' If non-`NULL`, an output file or connection to write the result to.
#'
#' @return [`StockholmMultipleAlignment`][inferrnal::StockholmMultipleAlignment-class]
#' object with modified RF line to mark the primer locations, or if `outfile` is
#' given, `NULL` invisibly.
#' @export
find_amplicon <- function(aln, fwd_primer, rev_primer,
                          trim = c("none", "retain", "remove"),
                          mark = TRUE, fwd_char = "{", 
                          rev_char = "}", outfile = NULL) {
  
  aln <- as_StockholmMSA(aln)
  
  if (is.character(fwd_primer)) fwd_primer <- Biostrings::DNAString(fwd_primer)
  assertthat::assert_that(methods::is(fwd_primer, "DNAString"))
  if (is.character(rev_primer)) rev_primer <- Biostrings::DNAString(rev_primer)
  assertthat::assert_that(methods::is(rev_primer, "DNAString"))
  rev_primer <- Biostrings::reverseComplement(rev_primer)
    
    assertthat::assert_that(
        assertthat::is.string(fwd_char),
        nchar(fwd_char) == 1L
    )
    
    assertthat::assert_that(
        assertthat::is.string(rev_char),
        nchar(rev_char) == 1L
    )
    
    assertthat::assert_that(fwd_char != rev_char)
  
  # find gaps in the reference alignment.  both "." and "-" are gaps
  gaps <- 
    mapply(
      Biostrings::union,
      Biostrings::vmatchPattern(".", aln@unmasked),
      Biostrings::vmatchPattern("-", aln@unmasked)
    )
  
  # get the ungapped reference sequences
  ungapped <-
    lapply(gaps, Biostrings::gaps, start = 1, end = ncol(aln)) |>
    mapply(
      FUN = function(x, ranges) x[ranges],
      x = aln@unmasked
    )
  if (methods::is(aln, "StockholmRNAMultipleAlignment")) {
    ungapped <- Biostrings::RNAStringSet(ungapped)
  }
  ungapped <- Biostrings::DNAStringSet(ungapped)
  
  # find the primers in the ungapped sequences, and map positions back into the
  # alignment
  result_fwd <- find_primer(ungapped, gaps, fwd_primer)
  result_rev <- find_primer(ungapped, gaps, rev_primer)
  
  if (isTRUE(mark)) {
      # replace the RF line with the new primer line.
      aln <- mark_ref_line(aln, result_fwd["start"], result_fwd["end"], fwd_char)
      aln <- mark_ref_line(aln, result_rev["start"], result_rev["end"], rev_char)
  }
  if (trim == "retain") {
      truncate_alignment(
          aln,
          outfile = outfile,
          start = result_fwd["start"],
          stop = result_rev["end"]
      )
  } else if (trim == "remove") {
      truncate_alignment(
          aln,
          outfile = outfile,
          start = result_fwd["end"] + 1L,
          stop = result_rev["start"] - 1L
      )
  } else if (!is.null(outfile)) {
      inferrnal::writeStockholmMultipleAlignment(aln, outfile)
      invisible(NULL)
  } else {
      aln
  }
}

modify_cm_rf <- function(infile, outfile, rf) {
  if (assertthat::is.readable(infile)) {
    infile <- file(infile, open = "rt")
  }
  assertthat::assert_that(
    methods::is(infile, "connection"),
    assertthat::is.string(rf)
  )
  if (assertthat::is.string(outfile)) {
    file.create(outfile)
    outfile <- file(outfile, open = "wt")
  }
  assertthat::assert_that(
    methods::is(outfile, "connection")
  )
  rf_width <- nchar(rf)
  
  while (length(l <- readLines(infile, 1000L)) > 0) {
    # make sure we have alignment mapping
    if (any(grepl("^MAP +no", l))) {
      stop("input CM does not have alignment mapping")
    }
    # if the CM didn't have an RF line before, it will when we're done with it.
    RF_lines <- which(grepl("^RF +no"))
    l[RF_lines] <- sub("no", "yes", l[RF_lines], fixed = TRUE)
    
    # modify RF characters for CM
    MAT_lines <- which(grepl("\\[ +MAT[PRL] ", l))
    for (i in MAT_lines) {
      fields <- strsplit(trimws(l[i]), " +")[[1]]
      spaces <- strsplit(l[i], "[^ ]+")[[1]]
      type <- fields[2]
      if (fields %in% c("MATL", "MATP")) {
        pos <- as.integer(fields[5])
        fields[9] <- substr(rf, pos, pos)
      }
      if (fields %in% c("MATR", "MATP")) {
        pos <- as.integer(fields[6])
        fields[10] <- substr(rf, pos, pos)
      }
      l[i] <- paste(c(spaces, fields)[order(c(seq_along(spaces), seq_along(fields)))], collapse = "")
    }
    
    # modify RF characters for HMM
    HMM_lines <- which(grepl(" *[1-9]\\d* +([0-9]\\.[0-9]+ +){4}[1-9][0-9]* +. +. +.", l))
    for (i in HMM_lines) {
      fields <- strsplit(trimws(l[i]), " +")[[1]]
      spaces <- strsplit(l[i], "[^ ]+")[[1]]
      pos <- fields[6]
      fields[8] <- substr(rf, pos, pos)
      l[i] <- paste(c(spaces, fields)[order(c(seq_along(spaces), seq_along(fields)))], collapse = "")
    }
    writeLines(outfile, l)
  }
  invisible(NULL)
}

extract_amplicon <- function(seqs, aln) {
  
}
