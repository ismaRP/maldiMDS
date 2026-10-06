

#' Estimate Q2E per peptide using a variant of a least squares regression
#'
#' @param peaks data.frame. Isotopic peaks data produced by [preprocess_spectra()]
#' @param lm_use Function. Model variant to use: [wlm_q2e_intercept()] (default), [lm_q2e_intercept],
#' [wlm_q2e()], or [lm_q2e()]
#' @param by_spectrum Logical. Whether to calculate a value per replicate/spectrum (`TRUE`, default)
#' or sample (`FALSE`).
#'
#' @returns A `data.frame` contaning Q2E estimates with the following columns:
#'  * `sample_name`
#'  * `spectrumId` optional, only if `by_spectrum` is `TRUE`
#'  * `pep_number`
#'  * `gamma_0` coefficient of the native, non-deamidated peptide
#'  * `gamma_1` coefficient of the peptide with 1 deamidation
#'  * `gamma_2` coefficient of the peptide with 2 demidations. Only for markers with
#'  2 glutamines, other with it is `NA`
#'  * `q2e` Ratio of deamidated to total. This is \deqn{1-\frac{gamma_0}{gamma_0 + gamma_1 + gamma_2}}
#'  * `intercept` optional, only if using a linear model with free intercept
#'  * `residual` Squared root of the sum of the squared residuals
#'  * `reliability` 1 - `residual`
#' @export
#' @importFrom dplyr group_by summarise mutate
wls_q2e = function(
    peaks,
    lm_use = c(wlm_q2e_intercept, lm_q2e_intercept, wlm_q2e, lm_q2e),
    by_spectrum = TRUE) {

  if (by_spectrum) {
    q2e_vals = peaks %>%
      group_by(sample_name, spectrumId, pep_number) %>%
      summarise(lm_use(norm_int, weight, deam_0, deam_1, deam_2)) %>%
      mutate(reliability = 1 - residual)
  } else {
    q2e_vals = peaks %>%
      group_by(sample_name, pep_number) %>%
      summarise(lm_use(norm_int, weight, deam_0, deam_1, deam_2)) %>%
      mutate(reliability = 1 - residual)
  }
  return(q2e_vals)

}



#' Calculate q2e for a single peptide and sample or replicate
#'
#' Fit a linear model from the theoretical native and deamidated envelopes to an
#' experimental isotopic envelope.
#'
#' For use after group_by on samples and peptide index, or providing the data
#' directly.
#'
#' @param norm_int Normalized intensity of the isotopic peaks
#' @param weight Weight of the peak as calculated by [preprocess_spectra].
#' If using a model without weighting, it is ignored
#' @param deam_0 Intensity of the peak in the native, non-deamidated theoretical envelope
#' @param deam_1 Itensity of the corresponding peak in the peptide with 1 deamidation
#' @param deam_2 Itensity of the corresponding peak in the peptide with 2 deamidations
#' (0 if the peptide only has one glutamine)
#' @param data If the previous are just column names in a data.frame, the data
#' must be passed here
#' @param return_model `logical`, whether to return the whole model
#'
#' @return a `data.frame` with `q2e`, intercept, gamma and residual estimates
#' @export
#' @importFrom stats lm na.exclude
#' @examples
#' \dontrun{
#' q2e_vals = peaks %>% filter(n_peaks > 0) %>%
#'   group_by(sample_name, spectrumId, pep_number) %>%
#'   summarise(wlm_q2e(norm_int, weight, deam_0, deam_1, deam_2))
#' }
lm_q2e_intercept = function(norm_int = NULL,
                            deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                            data = NULL,
                            return_model = FALSE) {
  if (!is.null(data)) {
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    intercept = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~deam_0+deam_1+deam_2,
      na.action = na.exclude, data = data)
    if (return_model) return(lm_model)
    intercept = lm_model$coefficients[1]
    gamma_0 = lm_model$coefficients[2]
    gamma_1 = lm_model$coefficients[3]
    gamma_2 = lm_model$coefficients[4]
    if (is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$residuals^2))/sum(!is.na(norm_int), na.rm = TRUE)
  }
  return(data.frame(q2e=q2e, intercept=intercept,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' @rdname lm_q2e_intercept
wlm_q2e_intercept = function(norm_int = NULL, weight = NULL,
                             deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                             data = NULL,
                             return_model = FALSE) {
  if (!is.null(data)) {
    data$weight = data$weight / sum(data$weight, na.rm = TRUE)
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    weight = weight /sum(weight, na.rm = TRUE)
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    intercept = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~deam_0+deam_1+deam_2,
      weights = weight, data=data,
      na.action = na.exclude)
    if (return_model) return(lm_model)
    intercept = lm_model$coefficients[1]
    gamma_0 = lm_model$coefficients[2]
    gamma_1 = lm_model$coefficients[3]
    gamma_2 = lm_model$coefficients[4]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$weights * lm_model$residuals^2))
  }
  return(data.frame(q2e=q2e, intercept=intercept,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' @rdname lm_q2e_intercept
wlm_q2e = function(norm_int = NULL, weight = NULL,
                   deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                   data = NULL,
                   return_model = FALSE) {
  if (!is.null(data)) {
    data$weight = data$weight / sum(data$weight, na.rm = TRUE)
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    weight = weight / sum(weight, na.rm = TRUE)
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~0+deam_0+deam_1+deam_2,
      weights = weight, data=data,
      na.action = na.exclude)
    if (return_model) return(lm_model)
    gamma_0 = lm_model$coefficients[1]
    gamma_1 = lm_model$coefficients[2]
    gamma_2 = lm_model$coefficients[3]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$weights * lm_model$residuals^2))
  }
  return(data.frame(q2e=q2e,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}


#' @rdname lm_q2e_intercept
lm_q2e = function(norm_int = NULL,
                  deam_0 = NULL, deam_1 = NULL, deam_2 = NULL,
                  data = NULL,
                  return_model = FALSE) {
  if (!is.null(data)) {
    if (all(is.na(data$norm_int))) allna = TRUE else allna = FALSE
  } else {
    if (all(is.na(norm_int))) allna = TRUE else allna = FALSE
  }
  if (allna) {
    q2e = NA
    gamma_0 = NA
    gamma_1 = NA
    gamma_2 = NA
    residual = NA
    if (return_model) return(NA)
  } else {
    lm_model = lm(
      norm_int~0+deam_0+deam_1+deam_2,
      na.action = na.exclude, data=data)
    if (return_model) return(lm_model)
    gamma_0 = lm_model$coefficients[1]
    gamma_1 = lm_model$coefficients[2]
    gamma_2 = lm_model$coefficients[3]
    if(is.na(gamma_2)) gamma_2 = 0
    q2e = 1 - (gamma_0/(gamma_0 + gamma_1 + gamma_2))
    residual = sqrt(sum(lm_model$residuals^2))/sum(!is.na(norm_int), na.rm = TRUE)
  }
  return(data.frame(q2e=q2e,
                    gamma_0=gamma_0, gamma_1=gamma_1, gamma_2=gamma_2,
                    residual=residual))
}



#' Estimation of peptide q2e in a single linear regression with nested coefficients
#'
#' Use a single call to "lm" with nested coefficients per spectra and peptide.
#' Effectively it will calculate separate coefficient for each spectra and peptide.
#'
#' WARNING: this cannot handle too large data.frames.
#' I will try to add `biglm` at some point.
#' Otherwise use group_by and the other functions
#'
#' @param peaks_df Data.frame with peaks data and theoretical isotopic envelopes
#' @param intensity Term used as the dependent variable. It could be the unscaled
#' intensity or the normalized at the spectra and peptide level
#' @param weights Column of weights in \code{peaks_df}
#' @param intercept Logical, whether to fit an intercept or not
#' @return A data.frame with q2e estimates
#' @importFrom stats lm as.formula
#' @export
#'
lm_q2e_oneshot = function(peaks_df, intensity='norm_int', weights=NULL, intercept=FALSE){

  # Prepare formula
  if (intercept) {
    lm_formula = sprintf(paste0(
      "%s ~ 0 + (spectra_name : pep_idx)",
        "+(spectra_name : pep_idx : deam_0)",
        "+(spectra_name : pep_idx : deam_1)",
        "+(spectra_name : pep_idx : deam_2)"),
      intensity)
  } else {
    lm_formula = sprintf(paste0(
      "%s ~ 0",
        "+(spectra_name : pep_idx : deam_0)",
        "+(spectra_name : pep_idx : deam_1)",
        "+(spectra_name : pep_idx : deam_2)"),
      intensity)
  }
  lm_model = lm(as.formula(lm_formula), data=peaks_df)

  result = lm_model$coefficients

  q2e_data = mapply(
    function(n, x){
      spectra_name = grep("spectra_name(.*)", n, value=T)
      spectra_name = strsplit(spectra_name, 'spectra_name')[[1]][2]
      pep_idx = grep("pep_idx\\d", n, value=T)
      pep_idx = strsplit(pep_idx, 'pep_idx')[[1]][2]
      ndeam = grep("deam_\\d", n, value=T)
      if (length(ndeam) == 0) ndeam = 'intercept'
      return(data.frame(spectra_name=spectra_name, pep_idx=pep_idx,
                        ndeam=ndeam, value=x))
    },
    strsplit(names(result), ":"),
    result, SIMPLIFY = F)

  q2e_data = do.call(rbind, q2e_data)
return(q2e_data)

}

