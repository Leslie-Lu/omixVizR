#' @title MR_2SLS
#' @description Perform Two-Stage Least Squares (2SLS) Mendelian Randomization analysis.
#' @param infile The input file containing the data.
#' @param outcome The outcome variable.
#' @param outcome_name The name of the outcome variable.
#' @param prs_cols_match The columns to match for PRS.
#' @param regexpr_pattern The regular expression pattern to use for matching.
#' @param standardise Whether to standardise the prs data, Default: TRUE
#' @param digits The number of decimal places to round to, Default: 3
#' @param .progress Whether to show progress, Default: TRUE
#' @return A list containing the results of the 2SLS analysis.
#' @details This function performs a Two-Stage Least Squares (2SLS) Mendelian Randomization analysis.
#' @seealso
#'  \code{\link[data.table]{fread}}, \code{\link[data.table]{as.data.table}}, \code{\link[data.table]{setattr}}, \code{\link[data.table]{setDT}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{reexports}}
#'  \code{\link[purrr]{map}}
#'  \code{\link[stringr]{str_extract}}
#'  \code{\link[Matrix]{nearPD}}
#'  \code{\link[corpcor]{cov.shrink}}
#' @rdname MR_2SLS
#' @export
#' @importFrom data.table fread as.data.table setnames setDT
#' @importFrom dplyr select contains
#' @importFrom purrr map_chr map map_lgl walk
#' @importFrom stringr str_extract
#' @importFrom Matrix nearPD
#' @importFrom corpcor cor.shrink
MR_2SLS = function(infile, outcome, outcome_name, prs_cols_match, regexpr_pattern, standardise = TRUE, digits = 3, .progress= TRUE) {
  message("Part1: Read and prepare the data...")
  data = if (is.character(infile)) {
    data.table::fread(infile, header = TRUE, stringsAsFactors = FALSE)
    } else {
    data.table::as.data.table(infile)
    }
  markers_prs= data %>% dplyr::select(dplyr::contains(prs_cols_match)) %>% colnames()
  markers = purrr::map_chr(markers_prs, ~ stringr::str_extract(.x, regexpr_pattern))
  markers_prs_v2 = gsub("-", "_", markers_prs, fixed = TRUE)
  markers_name = setNames(markers, markers_prs_v2)
  znorm= function(x){
    return(as.numeric(scale(x)))
  }
  ## --- Get parameters about the population ---
  count_valid= function(x){
    sum(is.finite(x))
  }
  effective_size= data[, sapply(.SD, count_valid), .SDcols=markers_prs]
  new_effective_size_names = gsub("-", "_", names(effective_size), fixed = TRUE)
  names(effective_size) = new_effective_size_names
  case_size    = data[ get(outcome) == 1, .N ]
  control_size = data[ get(outcome) == 0, .N ]

  ## --- Calculation of effective number of independent metabolites ---
  message("Part2: Calculate the effective number of indepedent metabolites...")
  resid_scale = function(biomarker){
    idx = grepl("-", biomarker, fixed = TRUE)
    other_cols= union(setdiff(colnames(data), c(markers_prs, outcome)), biomarker)
    covars_data= data[, .SD, .SDcols = other_cols]
    old_col_names = colnames(covars_data)
    new_col_names = gsub("-", "_", old_col_names, fixed = TRUE)
    data.table::setnames(covars_data, old_col_names, new_col_names)
    biomarker = gsub("-", "_", biomarker, fixed = TRUE)
    if (idx) {
      message("Processing biomarker (w/ -): ", markers_name[biomarker])
    }

    fml   = as.formula(sprintf("%s ~ .", biomarker))
    if (idx) {
      message("Model formula: ", deparse1(fml))
    }
    m = lm(fml, data = covars_data)
    z = resid(m)

    return(znorm(z))
  }
  resid_results = purrr::map(
    markers_prs,
    ~ resid_scale(.x),
    .progress = .progress
  ) %>% do.call(cbind, .) %>% as.data.frame()
  names(resid_results) = markers_prs
  message("The dimension of residuals matrix is: ", dim(resid_results)[1], " x ", dim(resid_results)[2])

  message("Checking the correlation matrix of residuals...")
  if (purrr::map_lgl(resid_results, ~ all(is.numeric(.x))) %>% all() == FALSE) {
    stop("Non-numeric values found in residuals matrix.")
  }
  nz_sd = purrr::map_lgl(resid_results, ~ sd(.x, na.rm = TRUE) > 0) # number of zero-sd columns
  enoughN_inf = purrr::map_lgl(resid_results, ~ count_valid(.x) > 3) # number of columns with enough non-infinite values
  keep = nz_sd & enoughN_inf
  if (sum(keep) >= 2) {
    message(paste0("Removing ", sum(!keep), " metabolites with zero standard deviation or insufficient non-missing/non-infinite values."))
    message("The removed metabolites are: ", paste(markers_prs[!keep], collapse = ", "))
    data.table::setDT(resid_results)
    resid_results[, (markers_prs[!keep]) := NULL]
    message("The new dimension of residuals matrix is: ", dim(resid_results)[1], " x ", dim(resid_results)[2])
  } else {
    stop("Not enough metabolites with non-zero standard deviation and sufficient non-missing/non-infinite values.")
  }
  message("Imputate missing values or infinite values with column medians...")
  resid_results_imputed = purrr::map(resid_results, ~ {
    x = .x
    x[is.infinite(x)] = NA_real_
    x[is.na(x)] = median(x, na.rm = TRUE)
    x = znorm(x)
    return(x)
  }) %>% do.call(cbind, .) %>% as.data.frame()
  message("The dimension of imputed residuals matrix is: ", dim(resid_results_imputed)[1], " x ", dim(resid_results_imputed)[2])

  corr_results = cor(resid_results_imputed, use = "everything", method = "pearson")
  corr_results = (corr_results + t(corr_results)) / 2
  corr_results = as.matrix(Matrix::nearPD(
    corr_results,
    corr      = TRUE,
    keepDiag  = TRUE,
    do2eigen  = TRUE
  )$mat)
  eig = tryCatch({
    eigen(corr_results, symmetric = TRUE, only.values = TRUE)$values
  }, error = function(e) {
    message("Error in eigen(corr_results): ", e$message,
            "  -> trying shrinkage correlation")
    corr_results = corpcor::cor.shrink(resid_results_imputed)
    corr_results = (corr_results + t(corr_results)) / 2
    corr_results = as.matrix(Matrix::nearPD(corr_results, corr=TRUE)$mat)
    eigen(corr_results, symmetric = TRUE, only.values = TRUE)$values
  })

  n_independent_metabolites = ((sum(eig)^2) / sum(eig^2)) %>% ceiling()
  message(paste0("The number of independent metabolites is ", n_independent_metabolites, " among ", length(markers_prs), " metabolites."))

  ## --- Stage 2 of Two stage least square (2SLS) estimation of MR ---
  message("Part3: Perform stage 2 of Two stage least square (2SLS) estimation of MR...")
  purrr::walk(markers_prs[!keep], ~ {
    message("The removed metabolite is: ", .x, " due to zero standard deviation or insufficient non-missing/non-infinite values.")
    print(table(data[[.x]], useNA = "always"))
  })
  data[, (markers_prs[!keep]) := NULL]
  message("The new dimension of data matrix is: ", dim(data)[1], " x ", dim(data)[2])
  markers_prs_keep = markers_prs[keep]
  if(standardise){
    data[, (markers_prs_keep) := lapply(.SD, znorm), .SDcols=markers_prs_keep]
  }
  univariate_mr_model = function(biomarker){
    idx = grepl("-", biomarker, fixed = TRUE)
    keep_cols = union(setdiff(colnames(data), c(markers_prs, biomarker)), biomarker)
    model_data = data[, .SD, .SDcols = keep_cols]
    old_col_names = colnames(model_data)
    new_col_names = gsub("-", "_", old_col_names, fixed = TRUE)
    data.table::setnames(model_data, old_col_names, new_col_names)

    x_cols = setdiff(colnames(model_data), outcome)
    biomarker = gsub("-", "_", biomarker, fixed = TRUE)
    if (idx) {
      message("Processing biomarker (w/ -): ", markers_name[biomarker])
    }
    fml   = as.formula(paste(outcome,"~", paste(x_cols, collapse = "+")))
    if (idx) {
      message("Model formula: ", deparse1(fml))
    }
    model = glm(fml, family= binomial(link = "logit"), data=model_data)

    return(model)
  }
  univariate_mr_results = purrr::map(
    markers_prs_keep,
    ~ lulab.utils::extract_logistic_model(
      univariate_mr_model(.x),
      markers_name,
      n_independent_metabolites,
      digits,
      effective_size,
      case_size,
      control_size,
      outcome_name
    ),
    .progress = .progress
  ) %>% do.call(rbind, .) %>% as.data.frame()
  message("The dimension of MR results is: ", dim(univariate_mr_results)[1], " x ", dim(univariate_mr_results)[2])

  out_results= list(
    correlation_matrix= corr_results,
    mr_results= univariate_mr_results
  )
  message("Done.")
  return(out_results)
}
