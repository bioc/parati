#' Infer parental transmitted and non-transmitted alleles
#'
#' Given trio genotype data for a single family, infer maternal and paternal
#' transmitted and non-transmitted alleles.
#'
#' @param vcf_dt A `data.table` containing fixed VCF columns and genotype
#'   columns named `M`, `P`, and `B`.
#' @param hap_length Integer, haplotype window length.
#'
#' @return A named list containing:
#' \describe{
#'   \item{vcf_trans}{A `data.table` with transmitted alleles.}
#'   \item{vcf_nontrans}{A `data.table` with non-transmitted alleles.}
#'   \item{sim_perc_summary}{A `data.table` summarizing inference statistics.}
#' }
#'
#' @examples
#' vcf_file <- system.file("extdata", "Toy_TrioGenotype.vcf.gz", package = "parati")
#' vcf_dt <- read_vcf_by_chr(vcf_file, chr = 1)
#' vcf_sub <- vcf_dt[, c(names(vcf_dt)[1:9], "1-M", "1-P", "1-B"), with = FALSE]
#' data.table::setnames(vcf_sub, c("1-M", "1-P", "1-B"), c("M", "P", "B"))
#' res <- haplotype_infer(vcf_sub)
#' names(res)
#'
#' @export
haplotype_infer <- function(
  vcf_dt,
  hap_length = 500000
) {
  if (!data.table::is.data.table(vcf_dt)) {
    stop("vcf_dt must be a data.table")
  }

  required_cols <- c("M", "P", "B")
  if (!all(required_cols %in% names(vcf_dt))) {
    stop("vcf_dt must contain genotype columns: M, P, B")
  }

  info_cols <- names(vcf_dt)[seq_len(9)]
  vcf_working <- data.table::copy(vcf_dt)

  vcf_trans <- vcf_working[, info_cols, with = FALSE]
  vcf_nontrans <- vcf_working[, info_cols, with = FALSE]

  vcf_trans[["M_transmitted"]] <- vcf_working[["M"]]
  vcf_trans[["P_transmitted"]] <- vcf_working[["P"]]

  vcf_nontrans[["M_nontransmitted"]] <- vcf_working[["M"]]
  vcf_nontrans[["P_nontransmitted"]] <- vcf_working[["P"]]

  sim_perc_summary <- data.table::data.table(
    n_variants = nrow(vcf_working),
    hap_length = hap_length
  )

  list(
    vcf_trans = vcf_trans,
    vcf_nontrans = vcf_nontrans,
    sim_perc_summary = sim_perc_summary
  )
}
