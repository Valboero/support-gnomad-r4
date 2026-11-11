#' @title Get variant population data
#' @inheritParams pkgParams
#' @param genomes The genome(s) associated with the given variant_id
#' @description
#' Returns the variant_id(s) for the alternate reference genome
#'
#' @return List of population data for the given variant_id(s)
#'
#' @examples
#' \dontrun{
#' getVariantPopData(varids = '1-88923261-A-C', genomes = "GRCh38")
#' }
#'
#' @import ghql
#' @importFrom glue glue glue_collapse
#' @importFrom jsonlite fromJSON
#' @export


getVariantPopData <- function(varids, genomes) {

  stopifnot(validGenomes(genomes))
  gmCon <- GraphqlClient$new(url = apiUrl())
  datasets <- getDatasets(genomes)

  qfmt <- '
    {genomes}_{gsub("-", "_", varids)}:
    variant(variantId: "{varids}", dataset: {datasets}) {{
      variant_id
      chrom
      pos

      exome {{
        populations {{
          id
          ac
          an
          homozygote_count
          hemizygote_count
        }}
      }}

      genome {{
        populations {{
          id
          ac
          an
          homozygote_count
          hemizygote_count
        }}
      }}
    }}
  '

  rsp <- .makeAndEvalQuery(qfmt, environment())

  if (is(rsp, "failedQuery")) return(rsp)

  resLst <- fromJSON(rsp, flatten = TRUE)$data

  procPopData <- function(x) {

    # Compatibility with different API field names
    varid <- x$variant_id

    if (is.null(varid)) {
      varid <- x$variantId
    }

    # Extract one sequencing type
    extractPop <- function(seq_type) {

      seq_data <- x[[seq_type]]

      if (is.null(seq_data) ||
          is.null(seq_data$populations)) {
        return(NULL)
      }

      pop <- as.data.frame(seq_data$populations)

      if (nrow(pop) == 0) return(NULL)

      out <- data.frame(
        varid = rep(varid, nrow(pop)),
        chrom = rep(x$chrom, nrow(pop)),
        pos = rep(x$pos, nrow(pop)),
        sequencing_type = seq_type,
        pop,
        row.names = NULL
      )

      # Calculate ancestry-specific allele frequency
      out$af <- ifelse(
        out$an > 0,
        out$ac / out$an,
        NA_real_
      )

      return(out)
    }

    # Keep exome and genome data separate
    results <- Filter(
      Negate(is.null),
      list(
        extractPop("exome"),
        extractPop("genome")
      )
    )

    if (length(results) == 0) return(NULL)

    do.call(rbind, results)
  }

  res <- lapply(resLst, procPopData)

  return(res)
}

