#' Example peptide table
#'
#' Peptide table from Nair et al. (2023)
#'
#' @format ## `parchment_peptides`
#' A data frame with 8 rows and 6 columns:
#' \describe{
#'   \item{pep_number}{Peptide number: pep1, pep2, ...}
#'   \item{pep_name}{
#'        Peptide name, in this case according to the ZooMS marker
#'        nomenclature}
#'   \item{mass}{M+H mass of the peptide}
#'   \item{sequence}{Peptide sequence}
#'   \item{n_hyp}{Number of hydroxyprolines}
#'   \item{label}{Optional additional label for the peptide, for example to show on plots}
#' }
#' @references
#' Nair B, Palomo IR, Markussen B, Wiuf C, Fiddyment S, Collins MJ (2023)
#' Parchment Glutamine Index (PQI): A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data.
#' Peer community journal, 3. https://doi.org/10.24072/pcjournal.230
"parchment_peptides"
