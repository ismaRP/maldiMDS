
#' Calculate MDS using linear mixed effect model from q2e estimates
#'
#' @param q2e_vals data.frame with q2e estimations per sample, replicate and peptide
#' @param logq Whether q values are log-scaled before entering the LME model
#' @param g Reliability power in the error normal distribution from linear mixed effects model
#' @param return_model `logical`. Whether to return the lme model object
#' @param outdir Optional. directory where results tables and plots are saved
#'
#' @return list contaning the following:
#' \itemize{
#'   \item \strong{sample}: data.frame with aggregated MDS estimates per sample
#'   \item \strong{pep}: data.frame of fitted and residuals per replicate and peptide
#'   \item \strong{estimates}: list of model estimates
#' }
#' \code{sample} data.frame has the following columns:
#' \itemize{
#'   \item Sample: sample name
#'   \item Prediction: sample random effects obtained from custom formula that. It is the log(MDS)
#'   \item sd: log(MDS) standard error
#'   \item RanefModel: sample random effects, as calculated by \code{\link[nlme]{ranef}}.
#'         It should be the same as Prediction, up to numerical precision
#'   \item MDS.Model, exponentials estimates. It is the MDS
#'
#' }
#'
#' \code{pep} data.frame contains:
#' \itemize{
#'   \item Sample, replicate and peptide IDs
#'   \item q: is the q calculated by the WLS from the isotopic distributions
#'   \item Reliability: minimum least square of q caluclated in the WLS step
#'   \item Loqq: log(q)
#'   \item resp: either q or log(q), according to whether logq was FALSE or TRUE
#'   \item Fitted: fitted q value
#'   \item Res: pearson residuals from Fitted
#'   \item Fitted0: fitted q values at the peptide level. It's equal to the peptide fixed effect
#'   \item Res0: pearson residuals from Fitted0
#' }
#'
#' \code{estimates} list contains
#' \itemize{
#'   \item alpha: peptide fixed effects
#'   \item sigma2_S: sample random effect variance
#'   \item sigma2_R: replicate random effect variance
#'   \item gamma: Reliability 2·gamma exponent
#'   \item sigma2: peptide fixed effects variance
#' }
#'
#' @importFrom nlme lme varComb varPower varIdent varFixed lmeControl
#' @importFrom nlme ranef
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom readr write_csv
#' @importFrom stats complete.cases fitted residuals
#' @export
#'
lme_mds = function(q2e_vals, logq=TRUE, g=NULL, outdir=NULL,
                   return_model=F){

  if (logq){
    q2e_vals = q2e_vals %>% mutate(resp = log(q2e))
  } else {
    q2e_vals = q2e_vals %>% mutate(resp = q2e)
  }

  if (is.null(g)) g = 'free'

  ## mixed effect model
  if (g == "free"){
    m = lme(
      resp~0+pep_number,
      random = ~1|sample_name/spectrumId,
      weights = varComb(
        varPower(-1/2, form = ~reliability),
        varIdent(form = ~1|pep_number)),
      control = lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data = q2e_vals)
  } else {
    m = lme(
      resp~0+pep_number,
      random = ~1|sample_name/spectrumId,
      weights = varComb(
        varFixed(~I(1/reliability)),
        varIdent(form = ~1|pep_number)
      ),
      control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
      data=q2e_vals)
    # m = lme(
    #   resp~0+Peptides,
    #   random=~1|Sample/Replicates,
    #   weights=varComb(varPower(fixed=g, form=~Reliability),
    #                   varIdent(form=~1|Peptides)),
    #   control=lmeControl(maxIter = 1000, msMaxIter = 1000, msMaxEval = 1000),
    #   data=q2e)
  }

  ## extract parameter estimates from lme object
  estimates_m = extract_estimates(m)
  if (g != "free"){
    estimates_m$gamma = g
  }
  prediction = predict_mds(m, estimates_m, logq=logq)
  ## mutate into the parent dataframe fitted values and residuals for plotting
  # q2e_m = q2e %>%
  # mutate(Fitted=fitted(m),
  #        Res=residuals(m, type = "pearson"),
  #        Fitted0=fitted(m, level = 0),
  #        Res0=residuals(m, level = 0, type = "pearson")) %>%
  # as_tibble()

  # q2e_m_pred = q2e_m %>% group_by(Sample) %>%
  #   summarise(predict_sample(
  #     Sample, Replicates, Peptides, Reliability, resp,
  #     estimates_m))

  ## here I have added untransformed and transformed prediction from both the model and function
  ## Prediction & MDS.PredictSample => from the function
  ## RanefModel & MDS.Model => from the model
  # if (logq){
  #   q2e_m_pred = q2e_m_pred %>%
  #     mutate(RanefModel=ranef(m)[['Sample']][,1],
  #            MDS.Model=exp(RanefModel),
  #            MDS.PredictSample=exp(Prediction),
  #            Sample = as.character(Sample))
  # } else {
  #   q2e_m_pred = q2e_m_pred %>%
  #     mutate(RanefModel=ranef(m)[['Sample']][,1],
  #            Sample = as.character(Sample),
  #            MDS.Model=RanefModel,
  #            MDS.PredictSample=Prediction)
  # }

  if (!is.null(outdir)){
    write_csv(
      prediction$pep,
      file.path(
        outdir,
        sprintf('MDS_pep_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
    write_csv(
      prediction$sample,
      file.path(
        outdir,
        sprintf('MDS_sample_estimates_gamma%3.3f.csv', estimates_m$gamma))
    )
  }

  prediction[['estimates']] = estimates_m
  if (return_model) prediction[['model']] = m

  return(prediction)
}


#' Predict MDS for samples
#'
#' @param model lme model object
#' @param new_q2e New q2e data for which MDS is predicted using the trained model
#' q is stored in a column called "resp"
#' @param logq whether q is in log-scale in q2e or new_q2e
#' @param estimates Contains estimates as extracted by [extract_estimates]:
#'
#' @importFrom nlme lme ranef
#' @importFrom stats predict
#' @importFrom dplyr mutate filter group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom stats complete.cases fitted residuals
#' @return `list` with peptide and sample level MDS predictions.
#' @export
predict_mds = function(model, estimates, new_q2e=NULL, logq=T){

  if (is.null(new_q2e)) {
    q2e = model$data
    mds = q2e %>% group_by(sample_name) %>%
      summarise(predict_sample(
        sample_name, spectrumId, pep_number, reliability, resp, estimates)) %>%
      mutate(MDS.Model = ranef(model)[['sample_name']][,1])
    if (logq) {
      mds = mds %>% mutate(
        # MDS.PredictSample = exp(MDS.PredictSample),
        MDS.Model = exp(MDS.Model),
        MDS.Manual = exp(MDS.Manual)
      )
    }
    q2e_m = q2e %>% ungroup() %>%
      mutate(Fitted = fitted(model),
             Res = residuals(model, type = "pearson"),
             Fitted0 = fitted(model, level = 0),
             Res0 = residuals(model, level = 0, type = "pearson")) %>%
      as_tibble()
    return(list('pep' = q2e_m, 'sample' = mds))
  } else {
    new_q2e = new_q2e %>%
      mutate(sample_name = as.factor(sample_name), spectrumId = as.factor(spectrumId),
             pep_number = as.factor(pep_number))
    # new_q2e = new_q2e %>% filter(residual>0) ## filter dataset to remove 0 reliabilities that resulted in Inf when taken the reciprocal

    mds = new_q2e %>% group_by(sample_name) %>%
      summarise(predict_sample(
        sample_name, spectrumId, pep_number, reliability, resp, estimates))
    q2e_m = new_q2e %>%
      mutate(
        predicted_q = predict(model, new_q2e)
      )
    if (logq){
      mds = mds %>% mutate(MDS.PredictSample = exp(MDS.PredictSample))
    }
    q2e_m = q2e_m %>%
      mutate(
        Pred = exp(predicted_q),
        Res = (resp - predicted_q)/sqrt(predicted_q))
    return(list('sample' = mds))
  }

}


#' Perform prediction of MDS and associated SD on a given sample
#'
#' @param sample_name Sample name
#' @param spectrumId Unique spectrum ID or replicate within a sample
#' @param pep_number Peptide identifier
#' @param reliability Reliability (1-residual)
#' @param resp Response variable, either q or Logq
#'
#' @param pars Contains the following estimates as extracted by [extract_estimates]:
#' \describe{
#'     \item{data_I_s}{part of the data frame corresponding to Sample s.}
#'     \item{alpha}{estimate for fixed effect}
#'     \item{sigma2_S}{estimate for variance for random effect of sample.}
#'     \item{sigma2_R}{estimate for variance for random effect of replication.}
#'     \item{gamma}{estimate for half power on Reliability}
#'     \item{sigma2}{estimate for variance for on the different peptides. This should be a named vector.}
#'}
#' @return tibble with predicted sample MDS and standard deviation
#' @export
#' @importFrom tibble tibble
predict_sample <- function(sample_name, spectrumId, pep_number, reliability, resp, pars) {

  data_I_s = data.frame(
    sample_name = sample_name, spectrumId = spectrumId, pep_number = pep_number,
    reliability = reliability, resp = resp)

  alpha = pars$alpha

  sigma2_S = pars$sigma2_S
  sigma2_R = pars$sigma2_R
  gamma = pars$gamma
  sigma2 = pars$sigma2

  if (nrow(data_I_s) == 0) {
    # prediction and variance without any data
    E_X = 0
    var_X = sigma2_S
  } else {
    # sanity check
    if (length(unique(data_I_s$sample)) != 1) stop("Observations must come from the same sample.")
    # build Xi
    Xi = sigma2_S * matrix(1, nrow(data_I_s), nrow(data_I_s)) +
      sigma2_R * outer(data_I_s$spectrumId, data_I_s$spectrumId, function(x,y){as.numeric(x == y)}) +
      diag((data_I_s$reliability^(2*gamma)) * sigma2[as.character(data_I_s$pep_number)], nrow = nrow(data_I_s))

    # Xi = sigma2_S * matrix(1, nrow(data_I_s), nrow(data_I_s)) +
    #   sigma2_R * outer(data_I_s$Replicates, data_I_s$Replicates, function(x,y){as.numeric(x == y)}) +
    #   diag((data_I_s$Reliability^(2*gamma)) * sigma2[as.character(data_I_s$Peptides)],nrow=nrow(data_I_s))
    # prediction
    # We don't need the estimate, as the ranef function provides it
    E_X = sigma2_S * sum(solve(Xi, data_I_s$resp - alpha[as.character(data_I_s$pep_number)]))
    # variance
    var_X = sigma2_S - (sigma2_S^2) * sum(solve(Xi, rep(1, nrow(data_I_s))))

  }

  # return
  names(var_X) = names(E_X) <- NULL
  # return(tibble(MDS.PredictSample = E_X, sd = sqrt(var_X)))
  return(tibble(MDS.Manual = E_X, sd = sqrt(var_X)))
}


#' Extract estimates from nlme model
#'
#' @param m nlme model
#'
#' @return list with estimates:
#' \describe{
#'     \item{data_I_s}{part of the data frame corresponding to Sample s.}
#'     \item{alpha}{estimate for fixed effect}
#'     \item{sigma2_S}{estimate for variance for random effect of sample.}
#'     \item{sigma2_R}{estimate for variance for random effect of replication.}
#'     \item{gamma}{estimate for half power on Reliability}
#'     \item{sigma2}{estimate for variance for on the different peptides. This should be a named vector.}
#'}
#' @importFrom nlme fixef
#' @importFrom stats coef
#'
extract_estimates = function(m) {
  q_data = m$data
  alpha.m = fixef(m)
  names(alpha.m) = substr(names(alpha.m), 11, nchar(names(alpha.m))) ## to align the code, remove Peptides from the names
  tmp = coef(m$modelStruct$reStruct, FALSE) * (m$sigma^2)
  sigma2_S.m = tmp["sample_name.var((Intercept))"] ## variance parameter estimate for sample
  sigma2_R.m = tmp["spectrumId.var((Intercept))"] ## variance parameter estimate for replicate
  gamma.m = coef(m$modelStruct$varStruct$A, FALSE)
  p1 = levels(q_data$pep_number)[1]
  sigma2.m = m$sigma * c(1, coef(m$modelStruct$varStruct$B, FALSE))^2
  names(sigma2.m) = c(p1, names(sigma2.m)[2:length(sigma2.m)])
  sigma2.m = sigma2.m[match(levels(q_data$pep_number), names(sigma2.m))]
  return(list(alpha = alpha.m, sigma2_S = sigma2_S.m, sigma2_R = sigma2_R.m,
              gamma = gamma.m, sigma2 = sigma2.m))
}



#' Quantile-Quantile plots of the linear mixed effect estimates per peptide
#'
#' @param mds_m data.frame of model estimates per replicate and peptide
#' @param title Plot title
#' @param pep_table A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column.
#' @param label_idx Index where to pull the labels from peptides
#' @param label_func labeller function to process labels. See \code{\link[ggplot2]{labeller}}
#' Default is label_value.
#'
#' @return `ggplot` plot object
#' @importFrom ggplot2 geom_qq geom_qq_line facet_wrap as_labeller labeller
#' @importFrom ggplot2 guide_legend guides aes
#' @importFrom ggplot2 ylab xlab ggtitle theme element_text theme_bw unit
#' @importFrom dplyr pull
#' @export
#'
pept_qqplot = function(mds_m, pep_table, title="", label_idx=2,
                       label_func = label_value){
  pept_labels = pull(pep_table, label_idx)
  pept_number = pull(pep_table, 1)
  names(pept_labels) = pept_number

  qq_plot = ggplot(mds_m) +
    geom_qq(aes(sample=Res, color = pep_number)) +
    geom_qq_line(aes(sample = Res, color = pep_number)) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~pep_number,
      labeller = labeller(pep_number = as_labeller(pept_labels, label_func))) +
    theme_bw() +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    ylab("quantiles of standardized residuals") +
    xlab("quantiles of standard normal") +
    ggtitle(title)
  return(qq_plot)
}

#' Fitted versus residuals plot per peptide
#'
#' @param mds_m data.frame of model estimates per replicate and peptide
#' @param title Plot title
#' @param pep_table A dataframe with peptide information. It must contain at least 3 columns,
#' peptide number or ID, name, and m/z. If NULL, default peptides are used.
#' The number or ID must have the form Pep# and be in the first column.
#' @param label_idx Index where to pull the labels from peptides
#' @param label_func labeller function to process labels. See \code{\link[ggplot2]{labeller}}
#' Default is label_value.
#'
#' @return `ggplot` plot object
#' @importFrom ggplot2 geom_point facet_wrap ggplot theme unit
#' @importFrom ggplot2 ylab xlab ggtitle aes element_text guides guide_legend
#' @importFrom dplyr pull
#' @export
#'
fvsr = function(mds_m, pep_table, title="", label_idx=2,
                label_func = label_value){

  pept_labels = pull(pep_table, label_idx)
  pept_number = pull(pep_table, 1)
  names(pept_labels) = pept_number

  fvsr_plot = ggplot(mds_m) +
    geom_point(aes(x = exp(Fitted), y = Res, color=pep_number),
               size = 2.5, alpha = 0.8) +
    # facet_wrap(~Peptides, scales="free") +
    facet_wrap(
      ~pep_number,
      labeller = labeller(pep_number = as_labeller(pept_labels, label_func))) +
    ylab("standardized residuals") +
    xlab("predicted q peptide") +
    theme_bw() +
    theme(legend.key.size=unit(1, "cm"),
          legend.text = element_text(size = 15),
          strip.text = element_text(size = 10)) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    ggtitle(title)
  return(fvsr_plot)
}



