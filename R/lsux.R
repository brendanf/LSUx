#' @importFrom magrittr "%>%" "%$%"
#' @importFrom dplyr .data
utils::globalVariables(c(".", "end", "start", "seq_from", "seq_to"))

#### LSUx ####
# functions to extract LSU and associated regions/domains from rRNA amplicons.
# also includes file I/O for stockholm and clustal alignments with RNA secondary
# structure information, as used by Infernal and Locarna.

#' Repair poor 5.8S detection in ITSx results using Infernal results.
#'
#' The 3' end of the 5.8S RNA is quite variable, and it is sometimes not
#' detected by the HMM used in ITSx.
#' Using a CMs, Infernal is able to more reliably delimit 5.8S.
#' This function uses the results of \code{\link[inferrnal]{cmsearch}} to
#' fill in missing 5.8S annotations in ITSx results.  It also updates the end of
#' ITS1 and beginning of ITS2 to match.
#'
#' @param itsx_result (\code{data.frame}) "positions" results from
#'        \code{\link[rITSx]{itsx}}
#' @param csearch_result (\code{data.frame}) results from
#'        \code{\link[inferrnal]{cmsearch}} run on the same set of sequences,
#'        using a 5.8S RNA model such as
#'        \href{https://rfam.xfam.org/family/RF00002/cm}{RF00002}
#'
#' @return a \code{\link[tibble]{tibble}} in the same format output by
#'   \code{\link[rITSx]{itsx}}, with updated positions for the boundaries of
#'   5.8S.
#' @export
#'
#' @examples
#' # load sample data from inferrnal
#' seqfile <- system.file(
#'     file.path("extdata", "sample.fasta"),
#'     package = "inferrnal"
#' )
#' # ITSx has trouble with some of the reads
#' seq <- Biostrings::readDNAStringSet(seqfile)[c(1,4,20,32,33,43,46,49)]
#' itsx_result <- rITSx::itsx(
#'     in_file = seq,
#'     positions = TRUE,
#'     complement = FALSE,
#'     cpu = 1,
#'     read_function = Biostrings::readDNAStringSet
#' )
#' pos <- itsx_result$positions
#' pos[pos$region == "5_8S",]
#' # find 5.8S using cmsearch
#' cm_5_8S <- system.file(
#'     file.path("extdata", "RF00002.cm"),
#'     package = "inferrnal"
#' )
#' cm_result <- inferrnal::cmsearch(seq, cm = cm_5_8S, cpu = 1)
#' # combine the results
#' merge_pos <- merge_5_8S(pos, cm_result)
#' merge_pos[merge_pos$region == "5_8S",]
merge_5_8S <- function(itsx_result, csearch_result) {
    csearch_result <- dplyr::select(
        csearch_result,
        seq_id = "target_name",
        "seq_from", "seq_to"
    ) %>%
        dplyr::mutate(region = "5_8S")

    # CMsearch is better at finding 5.8S than ITSx is at finding anything,
    # so csearch_result probably has hits in sequences that did not pass ITSx.
    # We don't have a comparable alternative to find LSU (it is so long that the
    # CM is much slower) so for now ignore those.
    out <- dplyr::full_join(itsx_result, csearch_result,
                            by = c("seq_id", "region"))

    # remove chimeras
    out <- out %>%
        dplyr::filter(
            is.na(comment) |
                stringr::str_detect(comment, "Chimer", negate = TRUE)
        ) %>%
        dplyr::group_by_at("seq_id") %>%
        dplyr::filter(sum(.data$region == "5_8S") <= 1) %>%
        dplyr::ungroup()

    # find the location difference between the ITSx HMM and the Rfam CM
    to_shift <- dplyr::filter(out, .data$region == "5_8S") %$%
        as.integer(stats::median(end - seq_to, na.rm = TRUE))
    from_shift <- dplyr::filter(out, .data$region == "5_8S") %$%
        as.integer(stats::median(start - seq_from, na.rm = TRUE))

    # When ITSx didn't find 5.8S, use the adjusted CM location
    out <-
        dplyr::mutate(
            out,
            start = dplyr::coalesce(
                .data$start,
                pmax(.data$seq_from + from_shift, 1L)
            ),
            end = dplyr::coalesce(
                .data$end,
                .data$seq_to + to_shift
            )
        )

    # Adjust ITS1 end and ITS2 start to match
    # This involves switching to "wide" format
    # Also remove comments about missing 5.8S if 5.8S has been found.
    out <- dplyr::full_join(
        dplyr::select(out, "seq_id", "length", "comment", "region", "start") %>%
            tidyr::spread(key = "region", value = "start"),
        dplyr::select(out, "seq_id", "length", "comment", "region", "end") %>%
            tidyr::spread(key = "region", value = "end"),
        by = c("seq_id", "length", "comment"),
        suffix = c("_start", "_end")
    ) %>%
        dplyr::mutate(
            ITS1_end = dplyr::coalesce(
                .data$`5_8S_start` - 1L,
                .data$ITS1_end
            ),
            ITS1_start = dplyr::coalesce(
                .data$ITS1_start,
                1L
            ),
            ITS2_start = dplyr::coalesce(
                .data$`5_8S_end` + 1L,
                .data$ITS2_start
            ),
            comment = ifelse(
                is.na(.data$`5_8S_start`) | is.na(.data$`5_8S_end`),
                .data$comment,
                stringr::str_replace(.data$comment, "[^!]*5.8S[^!]*! *", "")
            ) %>%
                dplyr::na_if("")
        )
    # Switch back to "long" format
    out <-
        dplyr::full_join(
            dplyr::select(out, "seq_id", "length", "comment",
                            tidyselect::ends_with("_start")) %>%
                dplyr::rename_all(stringr::str_replace, "_start", "") %>%
                tidyr::gather(key = "region", value = "start", -(seq.int(3))),
            dplyr::select(out, "seq_id", "length", "comment",
                            tidyselect::ends_with("_end")) %>%
                dplyr::rename_all(stringr::str_replace, "_end", "") %>%
                tidyr::gather(key = "region", value = "end", -(seq.int(3))),
            by = c("seq_id", "length", "comment", "region")
        )

    out
}

map_position <- function(alignment, x) {
    UseMethod("map_position")
}

map_position.MultipleAlignment <- function(alignment, x) {
    map_position.character(as.character(alignment), x)
}

map_position.character <- function(alignment, x) {
    if (length(x) == 1 && is.na(x)) return(rep(NA_integer_, length(alignment)))
    widths <- unique(nchar(alignment))
    assertthat::assert_that(
        assertthat::is.count(x),
        x >= 1,
        length(widths) == 1,
        x <= widths
    )
    if (x == 1L) return(rep(1L, length(alignment)))

    trimaln <- substr(alignment, 1L, x)
    gapcounts <- stringr::str_count(trimaln, "[.-]")
    pmax(x - gapcounts, 1L)
}

extract_rf_region <- function(rf, n, names) {
    assertthat::assert_that(
        assertthat::is.string(rf),
        is.character(n),
        all(nchar(n) == 1)
    )
    if (!missing(names)) {
        assertthat::assert_that(
            is.character(names),
            length(names) == length(n)
        )
    }
    out <- stringr::str_locate(rf, paste0(n, ".*", n))
    if (!missing(names)) {
        rownames(out) <- names
    }
    out
}

gap_free_width <- function(x, gapchars = ".-") {
    if (methods::is(x, "MultipleAlignment")) x <- x@unmasked
    if (is.character(x)) x <- Biostrings::BStringSet(x)
    assertthat::assert_that(methods::is(x, "XStringSet"))
    Biostrings::width(x) - c(Biostrings::letterFrequency(x, gapchars))
}

extract_LSU <- function(aln, rf, include_incomplete = FALSE, ...) {
    UseMethod("extract_LSU")
}

extract_LSU.MultipleAlignment <-
    function(aln, rf, include_incomplete = FALSE, ...) {

    extract_LSU.character(
        aln = as.character(aln@unmasked),
        rf = rf,
        include_incomplete = include_incomplete,
        seq_id = names(aln@unmasked),
        ...
    )
}

extract_LSU.character = function(aln, rf, include_incomplete = FALSE,
                                    seq_id = names(aln),
                                    length = gap_free_width(aln),
                                    ...) {
    limits <- extract_rf_region(
        rf, c(seq.int(9), LETTERS[seq.int(10)]),
        c("5_8S", paste0("LSU", seq.int(18)))
    )

    outhead <- tibble::tibble(
        seq_id = seq_id,
        length = length
    )

    out <- tibble::tibble(.rows = nrow(outhead))

    for (region in rownames(limits)) {
        startcol <- paste0(region, "_start")
        endcol <- paste0(region, "_end")

        out[[startcol]] <- map_position(aln, limits[region, "start"])
        out[[endcol]] <- map_position(aln, limits[region, "end"])

        out_of_range <- out[[startcol]] >= outhead$length &
            out[[endcol]] >= outhead$length
        out[[startcol]] <- ifelse(out_of_range, NA_integer_, out[[startcol]])
        out[[endcol]] <- ifelse(out_of_range, NA_integer_, out[[endcol]])
    }

    out[["ITS2_start"]] <- out[["5_8S_end"]] + 1L
    out[["ITS2_end"]] <- out[["LSU1_start"]] - 1L
    for (i in 2L:18L) {
        prename <- paste0("LSU", i - 1, "_end")
        postname <- paste0("LSU", i, "_start")
        if (isTRUE(include_incomplete) ||
            any(!is.na(out[[prename]]) & !is.na(out[[postname]]))) {
            zerosize <- out[[postname]] - out[[prename]] == 1L
            out[[paste0("V", i, "_start")]] <- ifelse(
                zerosize,
                NA_integer_,
                out[[prename]] + 1L
            )
            out[[paste0("V", i, "_end")]] <- ifelse(
                zerosize,
                NA_integer_,
                out[[postname]] - 1L
            )
        }
    }
    out <- purrr::discard(out, ~all(is.na(.)))
    out <- out[order(apply(out, 2, stats::median, na.rm = TRUE))]
    dplyr::bind_cols(outhead, out)
}

#' Extract variable/conserved domains from LSU rDNA
#'
#' Extracts alternating variable and conserved domains from the contiguous rDNA
#' regions which form the eukaryotic ribosomal large subunit,  i.e. 5.8S RNA,
#' ITS2, and 28S RNA. For the purposes of this document, this region will be
#' referred to as the 32S precursor RNA, as in humans, although its actual size
#' in Svedberg units varies between lineages.
#'
#' Input sequences should contain, at a minimum, a significant fraction of the
#' 5.8S RNA, which is used to define the 5' end of 32S. Any base pairs before
#' the 5' end of 5.8S will be considered to be ITS1 (\code{ITS1 = TRUE}) or
#' discarded (\code{ITS1 = FALSE}).  Input sequences should not extend past the
#' end of the 32S model at the 3' end.
#'
#' LSUx requires two covariance models: one for 5.8S, which is used in
#' \code{\link[inferrnal]{cmsearch}}, and one for 32S, which is used in
#' \code{\link[inferrnal]{cmalign}}.
#'
#' The 5.8S model can be
#' \href{https://rfam.xfam.org/family/RF00002}{RF00002 from Rfam} (the default),
#' or an equivalent.  It must be calibrated using \code{cmcalibrate} from
#' Infernal.
#'
#' The 32S model must include annotations in the reference line ("#=GC RF" in
#' the seed alignment) to distinguish conserved and variable regions. The
#' annotations should be sequential characters in the range \code{"1..9A..Z"}
#' for conserved domains, \code{"v"} for variable domains, and \code{"."} for
#' unaligned gaps in the seed alignment.
#' In the output, the conserved domains will be named "5_8S", "LSU1", "LSU2",
#' ...; the variable domains will be named "ITS2", "V1", "V2", ...
#'
#' Two example models are included,
#' both based on the
#' \href{https://github.com/rdpstaff/fungene_pipeline/blob/master/resources/RRNA_28S/model.cm}{RDP fungal LSU CM},
#' and annotated with variable regions according to Raué (1988).
#' The first, \code{system.file(file.path("extdata", "fungal_32S.cm"),
#' package = "LSUx")}, includes the full LSU region.  The second,
#' \code{system.file(file.path("extdata", "fungal_32S_LR5.cm"),
#' package = "LSUx")}, is truncated at the binding site of the LR5 primer, and
#' should be faster for input sequences which do not extend past that point.
#' The seed alignments are also provided.
#'
#' If generating similar truncated alignments with different endpoints, it is
#' critical to remove unpaired secondary structure elements from the
#' \code{"#=GC SS_cons"} line of the seed alignment.
#'
#' @param seq (filename readable by Infernal,
#'  \code{\link[Biostrings]{XStringSet}},
#'  \code{\link[ShortRead]{ShortRead}},
#'  or \code{character} vector) sequences
#'  to extract regions from
#' @param cm_5.8S (filename) covariance model for 5.8S rRNA
#' @param cm_32S (filename) covariance model for 32S pre-rRNA
#'      (5.8S, ITS2, and LSU)
#' @param glocal (\code{logical} scalar) if \code{TRUE}, use glocal alignment in
#'  \code{\link[inferrnal]{cmalign}}
#' @param ITS1 (\code{logical} scalar) if \code{TRUE}, include sequence fragment
#'  before 5.8S (if any) as ITS1
#' @param cpu (\code{integer} scalar) number of threads to use in Infernal calls
#'
#' @return a \code{\link[tibble]{tibble}} with one row for each region found for
#'  each input sequence.
#'  The columns are: \describe{
#'  \item{\code{seq_id} (\code{character})}{ the sequence name from
#'      \code{seq}}
#'  \item{\code{length} (\code{integer})}{ the length of the original sequence
#'       in base pairs}
#'  \item{\code{region} (\code{character})}{the name of the found domain}
#'  \item{\code{*_start} (\code{integer})}{the starting base for that domain in
#'    this sequence. Can be \code{"5_8S"}, \code{"ITS2"}, \code{"LSU1"},
#'    \code{"V2"}, \code{"LSU2"}, \code{"V3"}, etc.}
#'  \item{\code{end} (\code{integer})}{as \code{start}, but giving the end
#'      base for the domain.}}
#' @export
#'
#' @examples
#' # the sample data was amplified with primers ITS1 and LR5, so the truncated
#' # cm is appropriate.
#' seq <- system.file("extdata/sample.fasta", package = "inferrnal")
#' cm_32S_trunc <- system.file(
#'     file.path("extdata", "fungi_32S_LR5.cm"),
#'     package = "LSUx"
#' )
#' lsux(seq, cm_32S = cm_32S_trunc, ITS1 = TRUE, cpu = 1)
lsux <- function(
    seq,
    cm_5.8S = system.file(
        file.path("extdata", "RF00002.cm"),
        package = "inferrnal"
    ),
    cm_32S = system.file(
        file.path("extdata", "fungi_32S.cm"),
        package = "LSUx"
    ),
    glocal = TRUE,
    ITS1 = FALSE,
    cpu
) {
    assertthat::assert_that(
        assertthat::is.readable(cm_5.8S),
        assertthat::is.readable(cm_32S),
        assertthat::is.flag(glocal),
        assertthat::is.flag(ITS1)
    )

    seq <- protect_names(seq)

    futile.logger::flog.info("Beginning CM search.", name = "LSUx")
    cms <- inferrnal::cmsearch(
        cm = cm_5.8S,
        seq = seq$seq,
        glocal = glocal,
        cpu = cpu
    )

    # remove multiple hits
    cms <- dplyr::group_by_at(cms, "target_name")
    cms <- dplyr::filter(cms, all(.data$inc == "?") | .data$inc == "!")
    cms <- dplyr::filter(cms, dplyr::n() == 1)
    cms <- dplyr::ungroup(cms)
    futile.logger::flog.info(
        "%d sequences contained a single 5.8S hit.", nrow(cms),
        name = "LSUx"
    )
    if (is.character(seq)) {
        if (length(seq) == 1 && file.exists(seq)) {
            seq <- Biostrings::readBStringSet(seq)
        }
        else {
            seq <- Biostrings::BStringSet(seq)
        }
    }

    seq_idx <- as.integer(cms$target_name)
    seq_32S <- IRanges::narrow(seq$seq[seq_idx], start = cms$seq_from)

    futile.logger::flog.info("Beginning CM alignment.", name = "LSUx")
    aln <- inferrnal::cmalign(
        cmfile = cm_32S,
        seq = seq_32S,
        glocal = glocal,
        cpu = cpu
    )

    futile.logger::flog.info("Extracting LSU regions.", name = "LSUx")
    pos <- extract_LSU(aln = aln$alignment, rf = aln$RF)
    pos <- dplyr::mutate_at(
        pos,
        "seq_id",
        stringr::str_replace,
        "^[^|]*\\|",
        ""
    )
    pos <- dplyr::mutate_if(pos, is.integer, '+', cms$seq_from - 1L)
    if (isTRUE(ITS1)) {
        no_ITS1 <- pos$`5_8S_start` == 1 | is.na(pos$`5_8S_start`)
        if (all(no_ITS1)) {
            futile.logger::flog.warn(
                paste("ITS1 annotation was requested,",
                        "but no bases before 5.8S were found."),
                name = "LSUx"
            )
        } else {
            pos$ITS1_start <- ifelse(no_ITS1, NA_integer_, 1L)
            pos$ITS1_end <- ifelse(no_ITS1, NA_integer_, pos$`5_8S_start` - 1L)
            pos <- dplyr::select(
                pos,
                "seq_id",
                "length",
                "ITS1_start",
                "ITS1_end",
                dplyr::everything()
            )
        }
    }
    pos <- gather_regions(pos)
    pos$seq_id <- seq$seq_id[as.integer(pos$seq_id)]
    pos
}
