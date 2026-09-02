
#' Check whole number
#' 
#' @param x
#' @param tol numreric tolerance. Default is `.Machine$double.eps^0.5`
#' @return FALSE if `x` is not a whole number, TRUE if it is, given a numberic tolerance
is.wholenumber = function(x, tol = .Machine$double.eps^0.5) {
  if (!is.null(x)){
    abs(x - round(x)) < tol
  } else {
    FALSE
  }
}



#' Print only if verbose
#'
#' @param msg
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
print_progress = function(msg, verbose) {
  if (verbose) cat(msg)
}

#' Calculate the coefficients of each of the non or deamidated peptides
#' for a given global q2e
#'
#' @param q2e
#' @param max_Q
#'
#' @return
#' @export
#'
#' @examples
deamidation_coeffs = function(q2e, max_Q, ndeam) {
  triang = (max_Q+max_Q^2)/2
  unit_deam = (q2e)/triang
  coefs = lapply(max_Q, FUN=":",1)
  coefs = mapply(function(x, d) c((1-q2e), x*d), coefs, unit_deam)
  coefs = lapply(
    1:(ndeam+1),
    function(i) {
      sapply(coefs,
             function(x,i) {
               if(is.na(x[i])) 0 else x[i]
             }, i)
    })
  return(coefs)
}


#' Isotopic distribution of deamidated peptides
#' It produces the isotopic distribution of  peptides with a given extent
#' of deamidation.
#'
#' @param iso_peps
#' @param q2e
#'
#' @return
#' @export
#'
#' @examples
isotopic_deam_mat = function(iso_peps, max_Q, q2e=0.5, norm_func=NULL) {
  if (is.null(norm_func)) norm_func = function(x) return(1)
  ndeam = length(iso_peps) - 1
  coefs = deamidation_coeffs(q2e, max_Q, ndeam)
  abund = mapply(function(a,b) diag(a) %*% b, coefs, iso_peps, SIMPLIFY = F)
  iso_deam = Reduce('+', abund)
  iso_deam = iso_deam / apply(iso_deam, 1, max)
  # iso_peps[[paste0('deam', q2e)]] = iso_deam

  return(iso_deam)
}

#' Isotopic distribution of deamidated peptdes
#' It produces the isotopic distirbution of peptides with a given
#' extent of deamidation
#' @param iso_peps
#' @param q2e
#'
#' @return
#' @export
#'
#' @examples
isotopic_deam_df = function(iso_peps, q2e=0.5, norm_func=NULL, ...) {
  if (is.null(norm_func)) norm_func = function(x) return(1)
  iso_peps = iso_peps %>%
    group_by(pep_idx) %>%
    mutate(
      deam_comb = isotopic_deam(q2e=q2e, norm_func, max_Q, deam_0, deam_1, deam_2),
      q2e = q2e) %>%
    select(-c(deam_0, deam_1, deam_2))
  return(iso_peps)
}

#' Calculate the isotopic envelope of a deamidated peptide
#'
#' @param q2e Input
#' @param max_Q Number of deamidating Q sites
#' @param ... Relative intensity values of the isotopologues of the non-deamidated
#' versions of the peptide and the one with 1, 2,... deamidations
#'
#' @return
#' @export
#'
#' @examples
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
#' Calculates theoretical isotopic deistributions for a list of
#' sequences
#' @param seqs
#' @param nhyds Integer, vector
#' @param q2e Optional q2e value to calculate the theoretical distribution of
#' a peptide with this extent of deamidation
#'
#' @return Array of intensity distributions, without the masses
#' @export
#' @importFrom bacollite ms_iso
#' @importFrom stringr str_count
#'
#' @examples
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
