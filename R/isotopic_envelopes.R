
#' Calculates the chemical formula of a peptide sequence
#'
#' @param sequence `character`, peptide sequence
#'
#' @returns `character` with the formula using the Hill notation
#' @export
#'
#' @examples
#' peptide_formula(GVQGPPGPAGPR)
peptide_formula = function(sequence) {

  aa_formula = list(
    "A" = c(C=3, H=7,  N=1, O=2, S=0),
    "R" = c(C=6, H=14, N=4, O=2, S=0),
    "N" = c(C=4, H=8,  N=2, O=3, S=0),
    "D" = c(C=4, H=7,  N=1, O=4, S=0),
    "C" = c(C=3, H=7,  N=1, O=2, S=1),
    "E" = c(C=5, H=9,  N=1, O=4, S=0),
    "Q" = c(C=5, H=10, N=2, O=3, S=0),
    "G" = c(C=2, H=5,  N=1, O=2, S=0),
    "H" = c(C=6, H=9,  N=3, O=2, S=0),
    "I" = c(C=6, H=13, N=1, O=2, S=0),
    "L" = c(C=6, H=13, N=1, O=2, S=0),
    "K" = c(C=6, H=14, N=2, O=2, S=0),
    "M" = c(C=5, H=11, N=1, O=2, S=1),
    "F" = c(C=9, H=11, N=1, O=2, S=0),
    "P" = c(C=5, H=9,  N=1, O=2, S=0),
    "S" = c(C=3, H=7,  N=1, O=3, S=0),
    "T" = c(C=4, H=9,  N=1, O=3, S=0),
    "W" = c(C=11,H=12, N=2, O=2, S=0),
    "Y" = c(C=9, H=11, N=1, O=3, S=0),
    "V" = c(C=5, H=11, N=1, O=2, S=0)
  )

  aa = strsplit(toupper(sequence), "")[[1]]

  # Sum free amino-acid formulas
  chemform = Reduce(
    `+`,
    lapply(aa, function(x) aa_formula[[x]])
  )

  # Remove H2O for every peptide bond
  n_bonds = length(aa) - 1

  chemform["H"] = chemform["H"] - 2 * n_bonds
  chemform["O"] = chemform["O"] - n_bonds

  # Format in Hill notation
  chemform = chemform[chemform != 0]

  chemform = chemform[
    c(
      intersect(c("C", "H"), names(chemform)),
      sort(setdiff(names(chemform), c("C", "H")))
    )
  ]

  return(
    paste0(
      names(chemform),
      ifelse(chemform == 1, "", chemform),
      collapse = ""
      )
    )
}


#' Add modifications to a chemical formula, adding and substracting atoms
#'
#' @param pepform `character` initial chemical formula.
#' @param mods List of modifications to add to the peptide, with the number of
#' modifications, e.g. `list('deam'=1, 'hyd'=2)`
#'
#' @returns A `character` with the chemical formula in Hill notation of the
#' peptide with the modifications added
#' @export
#'
#' @examples
#' pepform = peptide_formula(GVQGPPGPAGPR)
#' modform = add_modifications(pepform, mods=list('deam'=1, 'hyd'=1))
add_modifications = function(pepform, mods=NULL) {
  if (is.null(mods)) return(pepform)
  for (m in names(mods)) {
    n = mods[[m]]
    if (n == 0) next
    if (!is.wholenumber(n)) stop('Non-integer number of modifications')
    # Add
    if (modifications[m,'formula_add'] != '') {
      pepform = MetaboCoreUtils::addElements(
        pepform,
        MetaboCoreUtils::multiplyElements(
          modifications[m,'formula_add'], n
        )
      )
    }
    # Subtract
    if (modifications[m,'formula_sub'] != '') {
      pepform = MetaboCoreUtils::subtractElements(
        pepform,
        MetaboCoreUtils::multiplyElements(
          modifications[m,'formula_sub'], n
        )
      )
    }
  }
  # Finally add [M+H]+ adduct
  pepform = MetaboCoreUtils::addElements(pepform, 'H')

  return(pepform)

}


#' Calculate number of isotopic substitutions / neutron additions
#'
#' @param iso_variants As retuened by [get_isotopic_variants]
#'
#' @returns An array with the number of netroon additions or isotopic substitutions
#' @export
#'
get_n_isosubs = function(iso_variants) {
  otherisotopes = c(
    '13C'=1,
    '2H'=1,
    '15N'=1,
    '18O'=2, '17O'=1,
    '34S'=2, '33S'=1, '36S'=4)
  isotopic_pos = apply(
    iso_variants, 1,
    function(iv) {
      sum(iv[names(otherisotopes)] * otherisotopes, na.rm = T)
    })
  return(isotopic_pos)
}


#' Wrapper around [enviPat::isopattern()]
#'
#' It is predefined to return probabilities and use `charge` 1.
#' @inheritParams enviPat::isopattern
#' @param pepform A single peptide checmical formula
#' @returns The same matrix as [enviPat::isopattern()]
#' @export
#'
#' @examples
#' pepform = peptide_formula(GVQGPPGPAGPR)
#' modform = add_modifications(pepform, mods=list('deam'=1, 'hyd'=1))
#' get_isotopic_variants(modform)
get_isotopic_variants = function(pepform, threshold=1e-4) {
  iso_variants = enviPat::isopattern(
    isotopes, pepform, charge = 1, verbose=FALSE, rel_to = 2, threshold=threshold)
  iso_variants = iso_variants[[pepform]]
  return(iso_variants)
}

#' Calculate the theoretical isotopic envelope
#'
#' @param n_isopeaks Number of isotopic peaks to calculate
#' @param threshold Probability threshold of the isotopic variants
#' @inheritParams peptide_formula
#' @inheritParams add_modifications
#' @returns A `matrix` with 3 columns: number of isotopic substitutions,
#' mass and probability. It has `n_isopeaks` columns.
#' @details
#' First it calculates all isotopic variants separately as [enviPat::isopattern()]
#' does. Then it calculates how many isotopic substitutions or neutrons are added
#' to each variant and aggregates isotopic variants with equivalent substitutions.
#' Only variants above the `threshold` are used. Decrease if needed to allow for
#' more `n_isopeaks`
#' @export
#'
#' @examples
#' iso_env = isotopic_envelope(
#'   'GPPGPQGPP', mods=list('deam'=1, 'hyd'=1), n_isopeaks = 6)
#'
isotopic_envelope = function(
    sequence,
    mods=NULL,
    n_isopeaks=5,
    threshold=1e-4
) {
  pepform = peptide_formula(sequence)
  pepform = add_modifications(pepform, mods)


  iso_variants = get_isotopic_variants(pepform, threshold=threshold)
  n_isosubs = get_n_isosubs(iso_variants)
  iso_variants = cbind(iso_variants, n_isosubs)
  colnames(iso_variants)[ncol(iso_variants)] = 'n_isosubs'

  agg_abundance = by(
    iso_variants[,'abundance'],
    iso_variants[,'n_isosubs'],
    sum, simplify = T)

  avg_mass = by(
    iso_variants[,c('m/z', 'abundance')],
    iso_variants[,'n_isosubs'],
    function(x) {
      wm = weighted.mean(x[,'m/z'], x[,'abundance'])
    }
  )

  if (length(avg_mass) < n_isopeaks) {
    stop('Not enough isotopic peaks, decrease threshold')
  }

  iso_env = cbind(
    iso_subs = unique(n_isosubs)[1:n_isopeaks],
    mass = avg_mass[1:n_isopeaks],
    prob = agg_abundance[1:n_isopeaks]
  )
  return(iso_env)
}

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
#' Calculates theoretical isotopic deistributions for a list of sequences.
#' @param seqs Array or list of peptide sequences
#' @param nhyds Integer, vector. Number of hydroxyprolines.
#' @param ndeam Maximum number of deamidated to consider. Usually the maximum number of Q's
#' in all sequences overall. It generates deamidated pepitde envelopes for this amount
#' of deamidations at most, even if a given peptide has more Q's.
#' @param n_isopeaks Number of isotopic peaks to generate
#' @param long_format `logival`, whether output a list of matrices with the isotopic envelopes
#' or a single `data.frame` (see Value section)
#' @param threshold See [enviPat::isopattern()].
#' @return If `long_format` is `TRUE`, it returns a long-format `data.frame` with
#' the relative intensity of each peptide and isotopic position for the native (non-deamidated),
#' and the verisions with 1 or more deamidations.
#' If `long_format` is `FALSE`, it returns a list with the isotopic patters for the native and
#' deamidated versions separated in each element of the list.
#' Either way, it returns probabilities as intensities.
#'
#' Array of intensity distributions, without the masses
#' @export
#' @importFrom bacollite ms_iso
#' @importFrom stringr str_count
#'
#' @examples
#' isodists = get_isodists(parchment_peptides$sequence, 3, parchment_peptides$n_hyp, sum, TRUE)
#' isotopic_deam_df(isodists, 0.7, sum)
get_isodists = function(seqs, nhyds, ndeam, n_isopeaks=5, long_format=F, threshold=1e-4){
  data("isotopes", package = "enviPat")
  if (!long_format) {
    iso_peps = list()
    for (d in 0:ndeam) {
      isotopes = array(0, dim = c(length(seqs), n_isopeaks))
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          mods = list('deam' = d, 'hyd' = nhyds[i])
          isodist = isotopic_envelope(seqs[i], mods=mods, n_isopeaks=n_isopeaks, threshold = threshold)
          isotopes[i, (d + 1):n_isopeaks] = isodist[1:(n_isopeaks - d), 3]
        }
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
      isotopes = numeric(nrow(iso_peps))
      for (i in 1:length(seqs)) {
        nQs = str_count(seqs[i], 'Q')
        if (nQs >= d) {
          mods = list('deam' = d, 'hyd' = nhyds[i])
          isodist = isotopic_envelope(seqs[i], mods=mods, n_isopeaks=n_isopeaks, threshold = threshold)
          isotopes[((i - 1)*n_isopeaks + 1 + d):(i*n_isopeaks)] = isodist[1:(n_isopeaks - d), 3]
        }
      }
      iso_peps[[paste0('deam_', d)]] = isotopes
    }
  }
  return(iso_peps)
}
