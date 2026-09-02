#' @title Preprocessing spectra for q2e estimation
#' @description
#' Performs smoothening, baseline removal and peak detection on MALDI samples.
#' From the peaks, isotopic peaks for a list of peptides are extracted.
#'
#' @details
#' Provide the input data either:
#' * `indir` path to spectra in mzML format
#' * provide a list of paths to `mzml_files`
#' * provide a `Spectra` object directly in `sps_mzr`.
#' If data is provided using more than one of the options, the `sps_mzr` is used, and then the `mzml_files`.
#' If metadata is provided, only the spectra specified in it will be used, even if `indir` contains more files
#' 
#' 
#' In order to split spectra by taxa to be analysed with different peptides:
#' * taxon_factor can be a factor specifying the split. The order must match that of the input data.
#'   Then use taxon_column to specify the column of the `peptides_user` that contains the taxon.
#' * If only `taxon_column` is specified, it must be a column in both metadata containing the taxonomic ID
#'   and in the `peptides_user` dataframe.
#'
#' @param indir Folder containing spectra in mzML format.
#' @param metadata Data frame with spectra metadata.
#' It contains the following columns:
#'  * `sample_name`
#'  * `replicate`
#'  *  Optionally, a column specifying taxonomic identification, as specified by `taxon_column`
#' Any spectra not present in `metadata` is not further analysed. Conversely, the function checks all the entries in `metadata`
#' have a corresponding file or spectrum.
#' @param mzml_files Paths to mzML files
#' @param spectrum_file_name If mzml_files are provided, whether to use file names
#' as spectra names. Otherwise, it is assumed the the spectra IDs are in the mzML
#' files' headers.
#' @param sps_mzr Spectra object
#' @param peptides_user DataFrame with peptides to be used to calculate glutamine q2e deamidation on.
#'        It must contain the following columns:
#'        * `mass` with peptide M+H monoisotopic masses
#'        * `sequence`
#'        * `n_hyp` number of hydroxyprolines
#'        * `pep_number` a peptide number or ID, which can be the same across taxa. Eg. pep1, pep2, ...
#'        * `pep_name` optional. Additional peptide name, this could be the ZooMS marker name
#'        * `taxon_column` optional. If For markers, taxa that this peptide marker belongs to, if the preprocessing
#'                         is split by taxa.
#' The default peptides are the ones from Nair et al. (2022).
#' The paper contains the details on the preprocessing procedure.
#' @param make_plot logical, whether to make a plot of the preprocessing step. Default is `FALSE`.
#'        If set to `TRUE` the input spectral data cannot contain more than 10 spectra. If you need to
#'        plot more spectra, call the function multiple times.
#' @param taxon_factor Factor to split spectra by taxon if a different set of peptides is
#'        to be used for each.
#' @param taxon_column Column in `metadata` and `peptides_user` with the taxonomic ID used
#'        for splitting the analysis.
#' @param verbose Whether to output progress
#' @param spectrumId_in_file logical. If spectra is given as one spectrum per mzML file and spectrumId is
#'        the file name, instead of being in the ID label of the mzML file. This is discouraged.
#' @param nreplicates integer, default NULL. Number of replicates per sample. If set, only samples with
#'                    `nreplicates` are analysed
#' @param ncores `integer`.
#' Number of cores used by the [Spectra::MsBackendMzR()] backend in [Spectra::peaksData()]
#' Default is `NULL`, in which case it uses `detectCores() - 2` with [BiocParallel::MulticoreParam()].
#' If set to 1, it uses [BiocParallel::SerialParam()].
#' @param chunk_size `integer`.
#' Processing chunk size i.e. how many spectra are loaded and processed in a chunk.
#' Keep in mind there will be \code{ncores} chunks processed in parallel.
#' Default is 40L.
#' @param ... MALDI-TOF Preprocessing and/or plotting parameters (see below).
#'    * `smooth_wma_hws` `integer`.
#'    Half-window size for WeightedMovingAverage smoothing method. Default 4.
#'    * `smooth_sg_hws` `integer`.
#'    Half-window size for SavitzkyGolay smoothing method. Default 6.
#'    * `iterations` `integer`.
#'    Iterations parameter for baseline detection using SNIP algorithm (see [MsCoreUtils::estimateBaseline()]). Default 50.
#'    * `halfWindowSize` `integer`.
#'    Half-window size parameter for local maximum detection. Default 20.
#'    * `snr` `numeric`.
#'    Signal-to-noise threshold above which peaks are considered.
#'    Only used for the SuperSmoother or MAD noise estimates.
#'    For local background estimation, see params `local_bg`, `mass_range` and `bg_cutoff` and `l_cutoff`.
#'    Default 2.
#'    * `k` `integer`.
#'    `k` parameter for [MsCoreUtils::refineCentroids()] Default is 0L, i.e. do not perform centroid refinement.
#'    * `threshold` `numeric`.
#'    threshold parameter for [MsCoreUtils::refineCentroids()] Default is 0.33, meaning only intensities above 0.33*maximal peak intenity are used.
#'    * `local_bg` `logical`.
#'    Whether to further to clean peaks of lists by modelling the local background.
#'    Default `FALSE`.
#'    See [MALDIzooMS::peaks_local_bg()].
#'    Ideally should work with a \code{snr} threshold of 0L.
#'    * `mass_range` `numeric`.
#'    Mass window to both sides of a peak to be considered for local background modelling
#'    Default 100.
#'    * `bg_cutoff` `numeric`.
#'    The peaks within the mass range with intensity below the \code{bg_cutoff} quantile
#'    are considered for background modelling. \code{bg_cutoff=1} keeps all peaks
#'    and \code{bg_cutoff=0.5} would only keep the bottom half. Default 0.5.
#'    * `l_cutoff` `numeric`.
#'    Likelihood threshold or p-value. Peaks with a probability of being modelled as
#'    background noise higher than this are filtered out.
#'    Default 1e-8.
#'    * `tolerance` `numeric`.
#'    Mass tolerance in Da between \code{mono_masses} with subsequent isotopic peaks
#'    and detected peaks. See [MsCoreUtils::closest()]
#'    Default 0.4.
#'    * `ppm` `numeric`.
#'    Parts-per-million added to tolerance. See [MsCoreUtils::closest()].
#'    Default 50.
#'    * `n_isopeaks` `integer`.
#'    Number of isotopic peaks to pick. Default is 5 and the maximum permitted.
#'    Default 5L.
#'    * `min_isopeaks` `integer`.
#'    If less than min_isopeaks consecutive (about 1 Da difference) isotopic peaks
#'    are detected, the whole isotopic envelope is discarded.
#'    Default 4L.
#'    * `norm_func` `function` or `NULL` (default).
#'    Function to normalize the isotopic distribution.
#'    If `NULL` (default), the isotopic peaks are normalised to the highest one.
#'    * `q2e` `numeric`. Only for plotting.
#'    A theoretical isotopic invelope with this `q2e` is overlaid in the plot in blue.
#'    * `peptide_labeller` `function`. Only for plotting.
#'    A function used to transform the facet labels. See [ggplot2::facet_wrap()] and [ggplot2::labeller()]
#' @return A list of dataframes, 1 per sample. Each dataframe has 3 columns,
#' m/z, intensity and signal-to-noise ratio for each of the n_isopeaks from each
#' peptide. Missing peaks are NAs.
#' @references
#' Nair, B. et al. (2022) ‘Parchment Glutamine Index (PQI):
#' A novel method to estimate glutamine deamidation levels in parchment collagen obtained from low-quality MALDI-TOF data’,
#' bioRxiv. doi:10.1101/2022.03.13.483627.
#' @importFrom MALDIzooMS get_spectra_name match_metadata_spectra separate_sample_replicate
#' @importFrom fs dir_ls
#' @importFrom dplyr mutate
#' @importFrom BiocParallel MulticoreParam SerialParam SnowParam register bpmapply
preprocess_spectra = function(
  indir=NULL,
  metadata=NULL,
  mzml_files=NULL,
  spectrum_file_name=NULL,
  sps_mzr=NULL,
  peptides_user=NULL,
  make_plot = FALSE,
  taxon_factor=NULL,
  taxon_column=NULL,
  verbose=FALSE,
  spectrumId_in_file=FALSE,
  nreplicates = NULL,
  ncores = NULL, chunk_size=40,
  ...) {

  # This function is the preprocessing interface

  if (all(is.null(indir), is.null(mzml_files), is.null(sps_mzr))){
    stop('Please provide data either as a Spectra object in sps_mzr, ',
         'a directory with mzML files in indir, ',
         'or a character vector with path(s) to mzML files in mzml_files')
  }

  if (is.null(sps_mzr)) {
    print_progress('Reading mzML headers into Spectra', verbose)
    if (is.null(mzml_files) & !is.null(indir)){
      mzml_files = dir_ls(indir)
    }
    sps_mzr = Spectra(mzml_files, source = MsBackendMzR(), centroided = FALSE)
    if (length(sps_mzr) == length(mzml_files) & spectrumId_in_file) {
      # Remove extension and use filename as spectrumId
      spectrum_id = path_ext_remove(path_file(sps_mzr$dataOrigin))
      sps_mzr$spectrumId = spectrum_id
    }
  }

  if (!is.null(metadata)) {
    nrows_metadata_old = nrow(metadata)
    nspectra_old = length(sps_mzr)
    tmp = match_metadata_spectra(metadata, sps_mzr)
    metadata = tmp[['metadata']]
    sps_mzr = tmp[['sps_mzr']]
    nrows_metadata_new = nrow(metadata)
    nspectra_new = length(sps_mzr)
    # TODO: if verbose write report on how many missing
    if (verbose) {
      md_loss = nrows_metadata_old - nrows_metadata_new
      sp_loss = nspectra_old - nspectra_new
      if (md_loss > 0 | sp_loss > 0){
        message(
          "\nThere are metadata entries without spectrum or spectra without metadata:\n",
          sprintf("Metadata entries without spectra: %d\n", md_loss),
          sprintf("Spectra without metadata entry: %d\n", sp_loss)
        )
      }
    }
  } else {
    metadata = separate_sample_replicate(sps_mzr$spectrumId)
    metadata$spectrum_id = sps_mzr$spectrumId
  }
  metadata = metadata %>% group_by(sample_name) %>%
    mutate(n_replicates = n()) %>% ungroup()

  if (is.wholenumber(nreplicates)) {
    keep = metadata$n_replicates >= nreplicates
    metadata = metadata[keep,]
    sps_mzr = sps_mzr[keep]
  }
  if (is.null(peptides_user)) {
    peptides_user = peptides
  }

  if (!is.null(taxon_column) & is.null(taxon_factor)) {
    if (is.null(metadata)) {
      stop(
        'Cannot split preprocessing by taxon.\n',
        sprintf('Please provide a metadata table with taxon_column %s', taxon_column))
    }
    taxon_factor = factor(metadata[[taxon_column]])
    peptide_factor = factor(peptides_user[[taxon_column]])
  } else if (!is.null(taxon_column) & !is.null(taxon_factor)) {
    taxon_factor = factor(taxon_factor)
    taxon_factor = droplevels(taxon_factor)
    peptide_factor = factor(peptides_user[[taxon_column]])
  } else {
    # Create a fake factor that won't split
    taxon_factor = factor(rep.int('unique_group', length(sps_mzr)))
    peptide_factor = factor(rep.int('unique_group', length(peptides_user)))
  }

  # Make sure levels of taxon_factor and peptide_factor match
  levels(taxon_factor) = sort(levels(taxon_factor))
  levels(peptide_factor) = sort(levels(peptide_factor))

  if (any(levels(taxon_factor) != levels(peptide_factor))) {    
    taxf_levels = paste0(levels(taxon_factor), collapse = '\n  -')
    taxf_levels = paste0('  -', taxf_levels)
  
    pepf_levels = paste0(levels(peptide_factor), collapse = '\n  -')
    pepf_levels = paste0('  -', pepf_levels)

    stop(
      'Levels of taxa for splitting peptides and spectral data ',
      'do not match\nor are in different order.\n',
      sprintf('Data levels:\n%s\nPeptide levels:\n%s\n', taxf_levels, pepf_levels)
    )
  }

  if (verbose) {
    if (length(levels(taxon_factor)) > 1) {
      levels_split = paste0(levels(taxon_factor), collapse = '\n  -')
      cat(
        sprintf('\nSplitting preprocessing by taxa:\n  -%s', levels_split)
      )
    } else {
      cat('Spectra not split by taxa.\n')
    }
  }

  if (make_plot) {
    print_progress('Setting ncores to 1 and processing in a single chunk.', verbose)
    ncores = 1
    chunk_size = length(sps_mzr)
  }

  if (is.null(ncores)) {
    ncores = detectCores() - 2
  } else if (ncores < 1) {
    ncores = detectCores() - 2
  }

  if (ncores == 1) {
    parparam = SerialParam(progressbar = verbose)
    print_progress('\nUsing 1 core in SerialParam', verbose)
  } else if (.Platform$OS.type == "windows") {
    parparam = SnowParam(workers=ncores, progressbar = verbose)
    print_progress(sprintf('\nUsing %s cores in SnowParam\n', ncores), verbose)
  } else {
    parparam = MulticoreParam(workers=ncores, progressbar = verbose)
    print_progress(sprintf('\nUsing %s cores in MulticoreParam\n', ncores), verbose)
  }
  
  if (!make_plot) {
    prep_by_group = function(sps_mzr_gr, peptides_gr, taxon) {
      print_progress(sprintf('\n\n______________\nProcessing %s\n', taxon), verbose)
      ## TODO: proprocessing function call
      peaks = .preprocess_spectra(sps_mzr = sps_mzr_gr, pep_table = peptides_gr,
        verbose = verbose, parparam=parparam, chunk_size=chunk_size, ...)
      print_progress('Done', verbose)

      peaks = peaks %>% mutate(taxon = taxon)

      return(peaks)
    }
    peaks_alltaxa = bpmapply(
      prep_by_group,
      split(sps_mzr, taxon_factor),
      split(peptides_user, peptide_factor),
      as.list(levels(taxon_factor)),
      SIMPLIFY = F,
      BPPARAM = SerialParam(progressbar = FALSE)
    )
    peaks_alltaxa = do.call(rbind, peaks_alltaxa)
  } else {
    .preprocessing_plot(sps_mzr = sps_mzr, pep_table = peptides_user, ...)
  }

}


#' Internal preprocessing function
#' 
#' Preprocess a Spectra object with MsBackendMzR
#' 
#' @importFrom MALDIzooMS smooth baseline_correction peak_detection
#' @importFrom MALDIzooMS peptide_pseudo_clusters peaks_local_bg
#' @importFrom parallel detectCores
#' @importFrom magrittr %>%
#' @importFrom Spectra Spectra addProcessing peaksData
#' @importFrom Spectra MsBackendMzR processingChunkSize
#' @importFrom BiocParallel MulticoreParam SerialParam SnowParam register
#' @importFrom dplyr rename bind_cols
.preprocess_spectra = function(
    sps_mzr,
    pep_table,
    verbose,
    parparam,
    chunk_size,
    smooth_wma_hws = 4L,
    smooth_sg_hws = 6L,
    iterations = 50L,
    halfWindowSize = 20L,
    snr = 2, k = 0L, threshold = 0.33,
    local_bg = FALSE,
    mass_range=100, bg_cutoff=0.5, l_cutoff=1e-8,
    tolerance = 0.4, ppm=50,
    n_isopeaks = 5,
    min_isopeaks = 4,
    norm_func = NULL){

  
  mono_masses = pep_table$mass
  
  register(parparam)

  processingChunkSize(sps_mzr) = chunk_size

  # Weighted Moving Average Smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'WeightedMovingAverage',
    hws = smooth_wma_hws, int_index = 'intensity', in_place = FALSE)
  # Savitzky-Golay Filter smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'SavitzkyGolay',
    hws = smooth_sg_hws, int_index = 'intensity', in_place = FALSE)

  
  # Baseline estimation on MA smoothed
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::baseline_correction, int_index = 'intensity_WeightedMovingAverage',
    keep_bl = FALSE, substract_index = 'intensity_SavitzkyGolay', in_place = TRUE,
    method = 'SNIP', iterations = iterations, decreasing = TRUE)
  # PEAKS
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::peak_detection, halfWindowSize = halfWindowSize,
    method = 'SuperSmoother', snr = snr, k = k, threshold = threshold,
    descending = TRUE, int_index = 'intensity_SavitzkyGolay',
    add_snr=TRUE)

  if (local_bg) {
    sps_mzr = addProcessing(
      sps_mzr, peaks_local_bg, mass_range = mass_range, bg_cutoff = bg_cutoff,
      l_cutoff = l_cutoff, int_index = 'intensity_SavitzkyGolay')
  }

  # GET ISOTOPIC CLUSTERS
  sps_mzr = addProcessing(
    sps_mzr, peptide_pseudo_clusters,
    mono_masses = mono_masses, n_isopeaks = n_isopeaks, min_isopeaks = min_isopeaks,
    tolerance = tolerance, ppm = ppm)


  print_progress('\tProcessing spectra ...\n', verbose)
  peaks = peaksData(sps_mzr, BPPARAM=param)
  names(peaks) = sps_mzr$spectrumId
  print_progress('\tPreparing peaks ... ', verbose)
  int_col = 'intensity_SavitzkyGolay'
  peaks = prepare_peaks(
    peaks, peptides_user = pep_table, n_isopeaks = n_isopeaks,
    int_column = int_col, norm_func=norm_func, q2e=NULL)
  print_progress('Done\n', verbose)
  return(peaks)

}


#' Preprocessing function for plotting
#' 
#' Preprocess a Spectra object so that it can be plotted
#' 
#' @param q2e `numeric`.
#' A theoretical isotopic invelope with this `q2e` is overlaid in the plot in blue.
#' @param peptide_labeller `function`.
#' A function used to transform the facet labels. See [ggplot2::facet_wrap()] and [ggplot2::labeller()]
#' @importFrom MALDIzooMS get_spectra_name separate_sample_replicate
#' @importFrom MALDIzooMS smooth baseline_correction peak_detection
#' @importFrom MALDIzooMS peptide_pseudo_clusters peaks_local_bg
#' @importFrom parallel detectCores
#' @importFrom magrittr %>%
#' @importFrom Spectra Spectra addProcessing peaksData
#' @importFrom Spectra MsBackendMzR processingChunkSize
#' @importFrom BiocParallel MulticoreParam SerialParam SnowParam register
#' @importFrom dplyr rename bind_cols
.preprocessing_plot = function(
    sps_mzr,
    pep_table,
    verbose,
    smooth_wma_hws = 4,
    smooth_sg_hws = 6,
    iterations = 50,
    halfWindowSize = 20,
    snr = 2, k = 0L, threshold = 0.33,
    local_bg = FALSE,
    mass_range=100, bg_cutoff=0.5, l_cutoff=1e-8,
    tolerance = 0.4, ppm=50,
    n_isopeaks = 5,
    min_isopeaks = 4,
    norm_func = NULL,
    q2e = NULL,
    peptide_labeller = NULL
  ) {

  register(SerialParam(progressbar = FALSE))

  # Weighted Moving Average Smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'WeightedMovingAverage',
    hws = smooth_wma_hws, int_index = 'intensity', in_place = FALSE)
  # Savitzky-Golay Filter smoothing
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::smooth, method = 'SavitzkyGolay',
    hws = smooth_sg_hws, int_index = 'intensity', in_place = FALSE)
  
  # Baseline estimation on MA smoothed
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::baseline_correction, int_index = 'intensity_WeightedMovingAverage',
    keep_bl = TRUE, substract_index = 'intensity_SavitzkyGolay', in_place = FALSE,
    method = 'SNIP', iterations = iterations, decreasing = TRUE)
  # Get spectra
  sp = peaksData(sps_mzr)
  names(sp) = sps_mzr$spectrumId
  sp = as.data.frame(sp) %>% rename(spectra_name = group_name)
  sp = bind_cols(sp, separate_sample_replicate(sp$spectra_name, sep = '_'))
  # PEAKS
  sps_mzr = addProcessing(
    sps_mzr, MALDIzooMS::peak_detection, halfWindowSize = halfWindowSize,
    method = 'SuperSmoother', snr = snr, k = k, threshold = threshold,
    descending = TRUE, int_index = 'intensity_SavitzkyGolay_bl_corr_SNIP',
    add_snr=TRUE)
  if (local_bg) {
    sps_mzr = addProcessing(
      sps_mzr, peaks_local_bg, mass_range = mass_range, bg_cutoff = bg_cutoff,
      l_cutoff = l_cutoff, int_index = 'intensity_SavitzkyGolay_bl_corr_SNIP')
  }

  # GET ISOTOPIC CLUSTERS
  sps_mzr = addProcessing(
    sps_mzr, peptide_pseudo_clusters,
    mono_masses = mono_masses, n_isopeaks = n_isopeaks, min_isopeaks = min_isopeaks,
    tolerance = tolerance, ppm = ppm)
  
  peaks = peaksData(sps_mzr, BPPARAM=param)
  names(peaks) = sps_mzr$spectrumId

  int_col = 'intensity_SavitzkyGolay_bl_corr_SNIP'
  peaks = prepare_peaks(
    peaks, peptides_user = peptides_user, n_isopeaks = n_isopeaks,
    int_column = int_col, norm_func=norm_func, q2e=q2e)
  
  .make_plot(
    sps_mzr, peaks, pep_table, n_isopeaks, norm_func = norm_func,
    peptide_labeller = peptide_labeller)
}


#' Normalize vector of intensities to the highest
#'
#' @param intensity Vector of intensities
#' @param norm_func Normalising function. Takes a vector of intensities as input.
#' @return
#' @export
#'
#' @examples
normalize_intensity = function(intensity, norm_func) {
  if (all(is.na(intensity))) {
    norm_int = rep(NA, length(intensity))
  } else {
    norm_int = intensity/norm_func(intensity[!is.na(intensity)])
  }
  return(norm_int)
}


#' Prepare list of peaks data into a data.frame
#'
#' @param peaks List of peaks matrix. Names are used as spectra name
#' @param n_isopeaks Number of isotopic peaks
#' @param peptides_user
#' A dataframe with peptide information. It must contain at least 3 columns,
#' \code{pep_number} or ID, \code{mass}, \code{sequence} and \code{h_hyp} (# of hydroxyprolines).
#' IF NULL, default are used, see details.
#' @param int_column Columns in peaks with the intensity to be used
#' @param norm_func Function to normalize the intensities of the isotopic envelope
#' @param q2e
#' If provided, it adds the theoretical isotopic distribution of peptides with
#' this extent of deamidation
#' @return A data.frame with isotopic peaks detected from data, theoretical isotopic
#' envelopes for 1 and 2 deamidations and other associated data.
#' @details The default peptides are the ones from Nair et al. (2022).
#' The paper contains the details on the preprocessing procedure.
#'
#' @importFrom dplyr mutate rename select ungroup arrange
#' @importFrom MALDIzooMS separate_sample_replicate
#' @export
#'
#' @examples
prepare_peaks = function(peaks, n_isopeaks, peptides_user, int_column='intensity',
                         norm_func=NULL, q2e=NULL) {


  if (is.null(norm_func)) norm_func = max
  peaks = as.data.frame(peaks) %>%
    rename(spectra_name = group_name)

  iso_peps = get_isodists(
    peptides_user$sequence, 2, peptides_user$n_hyp,
    norm_func=norm_func, long_format = T)

  if (!is.null(q2e)) {
    deam_iso = isotopic_deam_df(iso_peps, q2e, norm_func=norm_func)
    iso_peps = iso_peps %>%
      mutate(theor_deam = deam_iso$deam_comb)
  }

  n_spectra = length(unique(peaks$spectra_name))
  n_peptides = nrow(peptides_user)
  eps = 1e-5
  peaks = peaks %>%
    mutate(
      mass_pos = as.factor(rep(
        seq(n_isopeaks),
        times = n_spectra * n_peptides)),
      pep_idx = as.factor(rep(
        rep(seq(n_peptides), each = n_isopeaks),
        times = n_spectra)),
      pep_number = as.factor(rep(
        rep(peptides_user$pep_number, each = n_isopeaks),
        times = n_spectra)),
      ndeam = rep(
        rep(str_count(peptides_user$sequence, 'Q'), each = n_isopeaks),
        times = n_spectra),
      weight = SNR/sqrt(pmax(abs(delta_mass), eps))
      )



  if ('pep_name' %in% colnames(peptides_user)) {
    peaks = peaks %>%
      mutate(pep_name = as.factor(rep(
        rep(peptides_user$pep_name, each = n_isopeaks),
        times = n_spectra)))
  }

  # peaks = cbind(peaks, separate_sample_replicate(peaks$spectra_name))
  peaks = bind_cols(
    peaks, separate_sample_replicate(peaks$spectra_name, sep = '_'))

  peaks = Reduce(function(x, y) merge(x, y, all = TRUE, by = c('pep_idx', 'mass_pos')),
                 list(peaks, iso_peps))


  peaks[['intensity_use']] = peaks[[int_column]]
  peaks = peaks %>%
    arrange(sample_name, replicate, pep_idx, mass_pos) %>%
    group_by(spectra_name, pep_idx) %>%
    mutate(norm_int = normalize_intensity(intensity_use, norm_func=norm_func),
           n_peaks = sum(!is.na(intensity_use))) %>%
    ungroup() %>%
    # Transform NAs to 0?
    # mutate(intensity_use = replace(intensity_use, is.na(intensity_use), 0),
    #        norm_int = replace(norm_int, is.na(norm_int), 0)) %>%
    # Sanity check, make sure that if no peaks detected, all is NA
    mutate(intensity_use = replace(intensity_use, n_peaks == 0, NA),
           norm_int = replace(norm_int, n_peaks == 0, NA)) %>%
    select(-intensity_use)

  return(peaks)

}



#' Title
#'
#' @param x
#' @param n_isopeaks
#' @param min_isopeaks
#'
#' @return
#'
#' @examples
calc_n_frac_peaks = function(x, n_isopeaks, min_isopeaks) {
  fracs = list()
  for (n in c(min_isopeaks:n_isopeaks)) {
    fracs[[paste0('frac_', n)]] = sum(x == n)/length(x)
  }
  return(data.frame(fracs))
}

#' Plotting the number of isotopic peaks detected
#'
#' @param peaks
#' @param n_isopeaks
#' @param min_isopeaks
#' @param ...
#' @param marker_order
#'
#' @return
#' @export
#'
#' @examples
plot_n_peaks_per_peptide = function(peaks, n_isopeaks, min_isopeaks, marker_order, ...) {

  if (is.null(marker_order)) {
    marker_order = unique(peaks$pep_number)
  }

  a = peaks %>%
    group_by(pep_number, ...) %>%
    summarise(calc_n_frac_peaks(n_peaks, n_isopeaks, min_isopeaks)) %>%
    ungroup() %>%
    pivot_longer(cols = starts_with('frac'), names_to='n_of_peaks',
                 names_prefix = 'frac_', values_to='fraction') %>%
    mutate(
      n_of_peaks = factor(
        as.integer(n_of_peaks), levels=c(min_isopeaks:n_isopeaks)),
      pep_number = factor(pep_number, levels=marker_order))

  x = min_isopeaks:n_isopeaks
  mapped = (x - min(x)) / max(x - min(x)) * (9 - 1) + 1


  ggplot(a) +
    geom_col(aes(x=pep_number, y=fraction, fill=n_of_peaks),
             position='stack') +
    scale_fill_grey(
      '# of isotopic\npeaks',
      start=0.8, end=0.3) +
    ylab('Fraction of samples') +
    xlab('peptide') +
    scale_x_discrete(guide = guide_axis(n.dodge = 1, angle=45)) +
    facet_wrap(vars(...)) +
    theme_bw() +
    theme(
      panel.grid.minor=element_blank(),
      # axis.text.x=element_text(
      #   angle=40, vjust = 0.7, hjust=0.6, size=10))
      axis.text.x=element_text(size=10, hjust=1))
}






#' Make a faceted preprocessing plot per spectra and peptide
#'
#' @param sp
#' @param peaks
#' @param peptides_user
#' @param n_isopeaks
#'
#' @return
#' @export
#'
#' @examples
.make_plot = function(
    sp, peaks, peptides_user, n_isopeaks,
    norm_func=NULL, peptide_labeller) {

  if (is.null(norm_func)) norm_func = max
  peaks_mask = list()
  sp_mask = rep(NA, nrow(sp))
  for (i in seq_along(peptides_user$mass)) {
    mono_mz = peptides_user$mass[i]
    s = (sp$mz > (mono_mz-3)) & (sp$mz < (mono_mz+6))
    sp_mask[s] = peptides_user$pep_number[i]
    p = (peaks$mz > (mono_mz-3)) & (peaks$mz < (mono_mz+6))
    peaks_mask[[i]] = p
  }
  peaks_mask = Reduce('|', peaks_mask)
  peaks = peaks[peaks_mask,]

  sp$pep_number = sp_mask
  sp = sp[!is.na(sp$pep_number),]

  sort_idx = order(peptides_user$mass)
  sp = sp %>%
    group_by(spectra_name, pep_number) %>%
    mutate(norm_int = intensity/norm_func(intensity_SavitzkyGolay_bl_corr_SNIP, na.rm = TRUE),
           norm_int_wma = intensity_WeightedMovingAverage/norm_func(intensity_SavitzkyGolay_bl_corr_SNIP, na.rm = TRUE),
           norm_int_bl = baseline_SNIP/norm_func(intensity_SavitzkyGolay_bl_corr_SNIP, na.rm = TRUE),
           norm_int_bl_corr = intensity_SavitzkyGolay_bl_corr_SNIP/norm_func(intensity_SavitzkyGolay_bl_corr_SNIP, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(pep_number = factor(pep_number, levels=peptides_user$pep_number[sort_idx]))

  peaks = peaks %>%
    mutate(pep_number = factor(pep_number, levels=peptides_user$pep_number[sort_idx]))
  spp = ggplot(sp) +
    geom_vline(aes(xintercept=mz), data=peaks, color='grey40', linetype='dashed',
               linewidth=0.5) +
    # Raw int
    geom_line(aes(x = mz, y = norm_int), color='grey70', alpha=0.8, linewidth=0.8) +
    # Smooth int
    geom_line(aes(x = mz, y = norm_int_wma),
              color = 'grey10', alpha = 0.8) +
    # Baseline
    geom_line(aes(x = mz, y = norm_int_bl),
              color = 'blue1', linetype = "dashed", linewidth=0.8) +
    geom_line(aes(x = mz, y = norm_int_bl_corr),
              color = 'brown2', linetype = "solid", alpha=1, linewidth=1) +
    # geom_line(aes(x=mz, y=b_d), color='blue') +
    # geom_text(aes(label=QCflag), x=+Inf, y=+Inf, vjust=1.3, hjust=1.2,
    #           data=sps_clusters[sele]@backend@spectraData) +
    geom_point(aes(x = mz, y = norm_int), shape = 19, size = 2,
               data = peaks) +
    facet_grid(spectra_name~pep_number, scales = 'free', labeller=peptide_labeller) +
    ylab('Normalized intensity') +
    xlab('') +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size=12),
          panel.grid.major = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  if ('theor_deam' %in% colnames(peaks)) {
    spp = spp +
      geom_line(aes(x=mz, y=theor_deam), size=1, color='#26828EFF',
                          data = peaks)

  }
  return(spp)
}
