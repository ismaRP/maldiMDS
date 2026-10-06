#' Isotopic distribution of deamidated peptides
#' It produces the isotopic distribution of peptides with a given extent
#' of deamidation.
#'
#' @param iso_peps Isotopic envelopes
#' @param max_Q Maximum number of deamidations per peptide (equal to number of glutamines in the peptide)
#' @param q2e Fraction of deamidation
#' @param norm_func A function used to normalize isotopic envelope intensities
#' @param ... Vectors `"deam_0", "deam_1`, ... relative intensities of each ppeptide and isotopic
#' position for 0, 1, ... deamidations
#' @return
#' * `isotopic_deam_mat` returns `matrix` with a theoretical isotopic envelope
#' of a peptide per row, with the given q2e deamidation
#' * ``
#' @export
#'
#' @examples
#' isodists = get_isodists(parchment_peptides$sequence, 3, parchment_peptides$n_hyp, sum, TRUE)
#' isotopic_deam_df(isodists, 0.7, sum)
isotopic_deam_mat = function(iso_peps, max_Q, q2e=0.5, norm_func=NULL) {
  if (is.null(norm_func)) norm_func = function(x) return(1)
  ndeam = length(iso_peps) - 1
  coefs = deamidation_coeffs(q2e, max_Q, ndeam)
  abund = mapply(function(a,b) diag(a) %*% b, coefs, iso_peps, SIMPLIFY = F)
  iso_deam = Reduce('+', abund)
  iso_deam = iso_deam / apply(iso_deam, 1, norm_func)
  # iso_peps[[paste0('deam', q2e)]] = iso_deam

  return(iso_deam)
}

#' @rdname isotopic_deam_mat
#' @export
#' @importFrom dplyr pick starts_with
isotopic_deam_df = function(iso_peps, q2e=0.5, norm_func=NULL) {
  if (is.null(norm_func)) norm_func = function(x) return(1)
  iso_peps = iso_peps %>%
    group_by(pep_idx) %>%
    mutate(
      convolv_int = do.call(
        isotopic_deam,
        c(
          list(q2e, norm_func, max_Q),
          as.list(pick(starts_with('deam_')))
        )
      ),
      q2e = q2e) %>%
    select(-starts_with('deam'))
  return(iso_peps)
}


#' Calculate the coefficients of each of the non or deamidated peptides
#' for a given global q2e
#' @param ndeam Number of deamidation to consider
#' @describeIn isotopic_deam_mat Calculate the coefficients of each of the non or deamidated peptides
#' for a given global q2e
deamidation_coeffs = function(q2e, max_Q, ndeam) {
  triang = (max_Q+max_Q^2)/2
  unit_deam = (q2e)/triang
  coefs = lapply(max_Q, FUN=":",1)
  coefs = mapply(function(x, d) c((1-q2e), x*d), coefs, unit_deam)
  coefs = lapply(
    1:(ndeam+1),
    function(i) {
      sapply(coefs, function(x,i) if(is.na(x[i])) 0 else x[i], i)
    })
  return(coefs)
}

#' Calculate the isotopic envelope of a deamidated peptide
#'
#' @describeIn isotopic_deam_mat Calculate the isotopic envelope of a deamidated peptide
isotopic_deam = function(q2e=0.5, norm_func, max_Q, ...) {
  max_Q = max_Q[1]
  deam_mat = do.call(cbind, list(...))
  deam_mat = deam_mat[,1:(max_Q+1)]
  triang = (max_Q+max_Q^2)/2
  unit_deam = (q2e)/triang
  coefs = (max_Q:1)*unit_deam
  coefs = c((1-q2e), coefs)
  coefs = matrix(coefs, nrow = length(coefs))
  deam_comb = deam_mat %*% coefs
  deam_comb = deam_comb / norm_func(deam_comb)
  return(as.vector(deam_comb))
}


#' Get isotopic distributions
#' @description
#' Wrapper for bacollite ms_iso.
#' Calculates theoretical isotopic deistributions for a list of sequences.
#' @param seqs Array or list of peptide sequences
#' @param ndeam Maximum number of deamidated to consider. Usually the maximum number of Q's
#' in all sequences overall. It generates deamidated pepitde envelopes for this amount
#' of deamidations at most, even if a given peptide has more Q's.
#' @param nhyds Integer, vector. Number of hydroxyprolines.
#' @param norm_func Function to use for normalizing intensities, usually `sum` or `max`.
#' @param long_format `logival`, whether output a list of matrices with the isotopic envelopes
#' or a single `data.frame` (see Value section)
#' @return If `long_format` is `TRUE`, it returns a long-format `data.frame` with
#' the relative intensity of each peptide and isotopic position for the native (non-deamidated),
#' and the verisions with 1 or more deamidations.
#' If `long_format` is `FALSE`, it returns a list with the isotopic patters for the native and
#' deamidated versions separated in each element of the list.
#'
#' Array of intensity distributions, without the masses
#' @export
#' @importFrom bacollite ms_iso
#' @importFrom stringr str_count
#'
#' @examples
#' isodists = get_isodists(parchment_peptides$sequence, 3, parchment_peptides$n_hyp, sum, TRUE)
#' isotopic_deam_df(isodists, 0.7, sum)
get_isodists = function(seqs, ndeam, nhyds, norm_func, long_format=F){
  n_isopeaks = 5
  if (!long_format) {
    iso_peps = list()
    for (d in 0:ndeam) {
      isotopes = array(0, dim = c(length(seqs), n_isopeaks))
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          isodist = ms_iso(seqs[i], ndeamidations = d, nhydroxylations = nhyds[i])
          isodist[,2] = isodist[,2] / norm_func(isodist[,2])
          isotopes[i, (d + 1):n_isopeaks] = isodist[1:(n_isopeaks - d), 2]
        } # else {
        # isotopes[i,] = NaN
        # }
      }
      iso_peps[[paste0('deam_',d)]] = isotopes
    }
  } else {
    pep_idx = rep(1:length(seqs), each = n_isopeaks)
    mass_pos = rep(1:n_isopeaks, times = length(seqs))
    max_Q = rep(str_count(seqs, 'Q'), each = n_isopeaks)
    iso_peps = data.frame(
      pep_idx = pep_idx,
      mass_pos = mass_pos,
      max_Q = max_Q
    )
    for (d in 0:ndeam) {
      isotopes = array(0, dim = length(seqs)*n_isopeaks)
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          isodist = ms_iso(seqs[i], ndeamidations = d, nhydroxylations = nhyds[i])
          isodist[,2] = isodist[,2] / norm_func(isodist[,2])
          isotopes[((i - 1)*5 + 1 + d):(i*5)] = isodist[1:(n_isopeaks - d),2]
        } # else {
        # isotopes[((i-1)*5+1):(i*5)] = NaN
        # }
      }
      iso_peps[[paste0('deam_', d)]] = isotopes
    }
  }
  return(iso_peps)
}
