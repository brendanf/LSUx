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

#' @importClassesFrom Biostrings DNAStringSet
#' @importClassesFrom ShortRead ShortRead
methods::setAs(
    "DNAStringSet",
    "ShortRead",
    function(from) {
        ShortRead::ShortRead(
            magrittr::set_names(from, NULL),
            Biostrings::BStringSet(names(from))
        )
    }
)

methods::setAs(
    "ShortRead",
    "DNAStringSet",
    function(from) {
        out <- ShortRead::sread(from)
        names(out) <- as.character(ShortRead::id(from))
        out
    }
)

methods::setAs(
    "character",
    "ShortRead",
    function(from) {
        if (length(from) == 1 && file.exists(from)) {
            from <- tryCatch(
                ShortRead::readFasta(from),
                error = function(e) {
                    ShortRead::readFastq(from)
                }
            )
        } else {
            ShortRead::ShortRead(
                sread = Biostrings::DNAStringSet(from, use.names = FALSE),
                id = Biostrings::BStringSet(names(from))
            )
        }
    }
)

sreadq_to_qsDNAss <- function(from) {
    to = Biostrings::QualityScaledDNAStringSet(
        x = ShortRead::sread(from),
        quality = Biostrings::quality(from)
    )
    names(to) <- ShortRead::id(from)
    to
}


#' @importClassesFrom Biostrings QualityScaledXStringSet
#' @importClassesFrom ShortRead ShortReadQ
methods::setAs(
    "ShortReadQ",
    "QualityScaledXStringSet",
    sreadq_to_qsDNAss
)

methods::setAs(
    "ShortReadQ",
    "QualityScaledDNAStringSet",
    sreadq_to_qsDNAss
)

qsDNAss_to_sreadq <- function(seq) {
    ShortRead::ShortReadQ(
        sread = magrittr::set_names(methods::as(seq, "DNAStringSet"), NULL),
        quality = Biostrings::quality(seq),
        id = Biostrings::BStringSet(names(seq))
    )
}

#' @importClassesFrom Biostrings QualityScaledDNAStringSet
methods::setAs(
    "QualityScaledDNAStringSet",
    "ShortReadQ",
    qsDNAss_to_sreadq
)

methods::setAs(
    "QualityScaledDNAStringSet",
    "ShortRead",
    qsDNAss_to_sreadq
)

protect_names <- function(seq) {
  UseMethod("protect_names")
}

protect_names.ShortRead <- function(seq) {
  seq_id <- as.character(ShortRead::id(seq))
  seq@id <- Biostrings::BStringSet(as.character(seq_along(seq)))
  list(
    seq_id = seq_id,
    seq = seq
  )
}

protect_names.default <- function(seq) {
    seq_id <- names(seq)
    names(seq) <- as.character(seq_along(seq))
    list(
        seq_id = seq_id,
        seq = seq
    )
}

protect_names.character <- function(seq) {
    if (length(seq) == 1 && file.exists(seq) && endsWith(seq, ".fasta")) {
        seq <- Biostrings::readBStringSet(seq)
    } else {
        seq <- Biostrings::BStringSet(seq)
    }

    abc <- Biostrings::uniqueLetters(seq)
    if (all(abc %in% Biostrings::DNA_ALPHABET)) {
        seq <- Biostrings::DNAStringSet(seq)
        seq <- Biostrings::RNAStringSet(seq)
    } else if (all(abc %in% Biostrings::RNA_ALPHABET)) {
        seq <- Biostrings::RNAStringSet(seq)
    } else stop("Sequence alphabet should be DNA or RNA for LSUx.")
    protect_names.default(seq)
}
