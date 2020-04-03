gather_regions <- function(pos) {
  starts <- dplyr::select(pos, -tidyselect::ends_with("_end"))
  starts <- tidyr::gather(
    starts,
    key = "region",
    value = "start",
    tidyselect::ends_with("_start"),
    na.rm = TRUE
  )
  starts <- dplyr::mutate_at(
    starts,
    "region",
    stringr::str_replace,
    "_start$",
    ""
  )

  ends <- dplyr::select(pos, -tidyselect::ends_with("_start"))
  ends <- tidyr::gather(
    ends,
    key = "region",
    value = "end",
    tidyselect::ends_with("_end"),
    na.rm = TRUE
  )
  ends <- dplyr::mutate_at(ends, "region", stringr::str_replace, "_end$", "")

  hvars <- names(pos)
  hvars <- purrr::discard(hvars, endsWith, "_start")
  hvars <- purrr::discard(hvars, endsWith, "_end")
  joinvars <- c(hvars, "region")
  out <- dplyr::full_join(starts, ends, by = joinvars)
  dplyr::arrange(out, !!!rlang::parse_exprs(hvars), start)
}

spread_regions <- function(pos) {
  hvars <- setdiff(names(pos), c("region", "start", "end"))

  starts <- dplyr::select(pos, -"end")
  starts <- dplyr::mutate_at(starts, "region", paste0, "_start")
  starts <- tidyr::spread(starts, key = "region", value = "start")

  ends <- dplyr::select(pos, -start)
  ends <- dplyr::mutate_at(ends, "region", paste0, "_end")
  ends <- tidyr::spread(ends, key = "region", value = "end")

  out <- dplyr::full_join(starts, ends, by = hvars)
  outhead <- out[hvars]
  outvals <- dplyr::select(out, -!!hvars)
  outvals <- outvals[order(apply(outvals, 2, stats::median))]
  dplyr::bind_cols(outhead, outvals)
}

as_sread <- function(seq) {
  UseMethod("as_sread")
}

as_sread.ShortRead <- function(seq) {
  seq
}

as_sread.DNAStringSet <- function(seq) {
  ShortRead::ShortRead(
    magrittr::set_names(seq, NULL),
    Biostrings::BStringSet(names(seq))
  )
}

as_sread.character <- function(seq) {
  if (length(seq) == 1 && file.exists(seq)) {
    seq <- tryCatch(
      ShortRead::readFasta(seq),
      error = function(e) {
        ShortRead::readFastq(seq)
      }
    )
  } else {
    ShortRead::ShortRead(
      sread = Biostrings::DNAStringSet(seq, use.names = FALSE),
      id = Biostrings::BStringSet(names(seq))
    )
  }
}

as_sread.QualityScaledDNAStringSet <- function(seq) {
  ShortRead::ShortReadQ(
    sread = magrittr::set_names(methods::as(seq, "DNAStringSet"), NULL),
    quality = Biostrings::quality(seq),
    id = Biostrings::BStringSet(names(seq))
  )
}

name_protected_sread <- function(seq) {
  UseMethod("name_protected_sread")
}

name_protected_sread.ShortRead <- function(seq) {
  seq_id <- as.character(ShortRead::id(seq))
  seq@id <- Biostrings::BStringSet(as.character(seq_along(seq)))
  list(
    seq_id = seq_id,
    seq = seq
  )
}

name_protected_sread.default <- function(seq) {
  name_protected_sread.ShortRead(as_sread(seq))
}
