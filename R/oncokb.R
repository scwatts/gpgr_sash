#' @export
read_oncokb <- function(x) {
  readr::read_tsv(x) |>
    dplyr::filter(
      `OncoKB Annotated` == "Yes"
    ) |>
    dplyr::pull("Hugo Symbol")
}

#' Unsurprisingly this is inefficient in R, we get ~4,000 rows/min
#' @export
set_gene_oncokb_links <- function(x) {

  delimiters <- ",&-"
  delimiter_re <- paste0("[", delimiters, "]")

  p <- stringr::str_subset(oncokb_genes, "-") |>
    purrr::map(function(n) {
      r <- paste0("(?:^|", delimiter_re, ")", n, "(?:", delimiter_re, "|$)")
      stringr::str_locate_all(x, r) |> purrr::pluck(1) |> tibble::as_tibble()
    }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(start, end) |>
    dplyr::mutate(
      end = ifelse(end == nchar(x), end, end - 1),
      processed = TRUE,
    )

  if (nrow(p) < 1) {
    d <- tibble::tibble(start=1, end=nchar(x), processed=FALSE)
  } else {

    d <- p
    last <- 0
    n <- 0

    for (i in 1:nrow(p)) {
      if ((last + 1) != p$start[i]) {
        d <- dplyr::add_row(d, start=last+1, end=p$start[i]-1, processed=FALSE, .before=i+n)
        n <- n + 1
      }
      last <- p$end[i]
    }

    if (p$end[i] != nchar(x)) {
      d <- dplyr::add_row(d, start=last+1, end=nchar(x), processed=FALSE)
    }
  }

  d |>
    dplyr::mutate(
      token = stringr::str_sub(x, start, end),
    ) |>
    dplyr::select(c(token, processed)) |>
    purrr::pmap(function(token, processed) {

      if (!processed) {
        tokens <- stringr::str_split(token, paste0("(?=", delimiter_re, ")"), simplify=TRUE) |>
          purrr::discard(~ nchar(.x) == 0 )
      } else {
        tokens <- list(token)
      }

      tokens |> purrr::map(function(t) {

        ld <- NULL
        ln <- NULL

        if (stringr::str_detect(t, paste0("^", delimiter_re))) {
          ld <- stringr::str_sub(t, 1, 1)
          ln <- stringr::str_sub(t, 2, nchar(t))
        } else {
          ln <- t
        }

        if (ln %in% oncokb_genes) {
          ln <- paste0("<a href='https://www.oncokb.org/gene/", ln, "' target='_blank'>", ln, "</a>")
        }

        return(c(ld, ln))

      })
    }) |>
    unlist() |>
    stringr::str_c(collapse="")
}
