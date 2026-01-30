# =============================================================================
# FIXED: R/forest_search_revised.R - forestsearch() roxygen2 documentation
# =============================================================================
# 
# This file contains the COMPLETE corrected roxygen2 documentation block for
# the forestsearch() function. The fix adds the missing @param stop_threshold.
#
# INSTRUCTIONS:
# 1. Open R/forest_search_revised.R
# 2. Find the roxygen2 documentation block starting with:
#    #' @title ForestSearch: Exploratory Subgroup Identification
# 3. Replace the entire @param section with the content below
#
# =============================================================================

#' @title ForestSearch: Exploratory Subgroup Identification
#'
#' @description
#' Identifies subgroups with differential treatment effects in clinical trials
#' using a combination of Generalized Random Forests (GRF), LASSO variable
#' selection, and exhaustive combinatorial search with split-sample validation.
#'
#' @param df.analysis Data frame. Analysis dataset with required columns.
#' @param outcome.name Character. Name of time-to-event outcome variable. Default "tte".
#' @param event.name Character. Name of event indicator (1=event, 0=censored). Default "event".
#' @param treat.name Character. Name of treatment variable (1=treatment, 0=control). Default "treat".
#' @param id.name Character. Name of subject ID variable. Default "id".
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param flag_harm.name Character. Name of true harm flag for simulations (optional).
#' @param confounders.name Character vector. Names of candidate subgroup-defining variables.
#' @param parallel_args List. Parallel processing configuration:
#'   \describe{
#'     \item{plan}{Character. One of "multisession", "multicore", "callr", "sequential"}
#'     \item{workers}{Integer. Number of parallel workers}
#'     \item{show_message}{Logical. Show parallel setup messages}
#'   }
#' @param df.predict Data frame. Prediction dataset (optional).
#' @param df.test Data frame. Test dataset (optional).
#' @param is.RCT Logical. Is this a randomized controlled trial? Default TRUE.
#' @param seedit Integer. Random seed. Default 8316951.
#' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".
#' @param use_lasso Logical. Use LASSO for variable selection. Default TRUE.
#' @param use_grf Logical. Use GRF for variable importance. Default TRUE.
#' @param grf_res GRF results object (optional, for reuse).
#' @param grf_cuts List. Custom GRF cut points (optional).
#' @param max_n_confounders Integer. Maximum confounders to consider. Default 1000.
#' @param grf_depth Integer. GRF tree depth. Default 2.
#' @param dmin.grf Integer. Minimum events for GRF. Default 12.
#' @param frac.tau Numeric. Fraction of tau for RMST. Default 0.6.
#' @param conf_force Character vector. Variables to force include (optional).
#' @param defaultcut_names Character vector. Default cut variable names (optional).
#' @param cut_type Character. Cut type ("default" or "custom"). Default "default".
#' @param exclude_cuts Character vector. Variables to exclude from cutting (optional).
#' @param replace_med_grf Logical. Replace median with GRF cuts. Default FALSE.
#' @param cont.cutoff Integer. Cutoff for continuous vs categorical. Default 4.
#' @param conf.cont_medians Named numeric vector. Median values for continuous variables (optional).
#' @param conf.cont_medians_force Named numeric vector. Forced median values (optional).
#' @param n.min Integer. Minimum subgroup size. Default 60.
#' @param hr.threshold Numeric. Minimum HR for candidate subgroups. Default 1.25.
#' @param hr.consistency Numeric. Minimum HR for consistency validation. Default 1.0.
#' @param sg_focus Character. Subgroup selection focus. One of "hr", "hrMaxSG", "maxSG",
#'   "hrMinSG", "minSG". Default "hr".
#' @param fs.splits Integer. Number of splits for consistency evaluation (or maximum
#'   splits when \code{use_twostage = TRUE}). Default 1000.
#' @param m1.threshold Numeric. Maximum median survival threshold. Default Inf.
#' @param pconsistency.threshold Numeric. Minimum consistency proportion. Default 0.90.
#' @param stop_threshold Numeric. Early stopping threshold for consistency
#'   evaluation. When a candidate subgroup's estimated consistency probability
#'   exceeds this threshold, evaluation stops early. Default 0.95.
#' @param showten_subgroups Logical. Show top 10 subgroups. Default FALSE.
#' @param d0.min Integer. Minimum control arm events. Default 12.
#' @param d1.min Integer. Minimum treatment arm events. Default 12.
#' @param max.minutes Numeric. Maximum search time in minutes. Default 3.
#' @param minp Numeric. Minimum prevalence threshold. Default 0.025.
#' @param details Logical. Print progress details. Default FALSE.
#' @param maxk Integer. Maximum number of factors per subgroup. Default 2.
#' @param by.risk Integer. Risk table interval. Default 12.
#' @param plot.sg Logical. Plot subgroup survival curves. Default FALSE.
#' @param plot.grf Logical. Plot GRF results. Default FALSE.
#' @param max_subgroups_search Integer. Maximum subgroups to evaluate. Default 10.
#' @param vi.grf.min Numeric. Minimum GRF variable importance. Default -0.2.
#' @param use_twostage Logical. Use two-stage sequential consistency algorithm for
#'   improved performance. Default FALSE for backward compatibility. When TRUE,
#'   \code{fs.splits} becomes the maximum number of splits, and early stopping
#'   is enabled. See Details.
#' @param twostage_args List. Parameters for two-stage algorithm (only used when
#'   \code{use_twostage = TRUE}):
#'   \describe{
#'     \item{n.splits.screen}{Integer. Splits for Stage 1 screening. Default 30.}
#'     \item{screen.threshold}{Numeric. Consistency threshold for Stage 1. Default
#'       is automatically calculated to provide ~2.5 SE margin.}
#'     \item{batch.size}{Integer. Splits per batch in Stage 2. Default 20.}
#'     \item{conf.level}{Numeric. Confidence level for early stopping. Default 0.95.}
#'     \item{min.valid.screen}{Integer. Minimum valid Stage 1 splits. Default 10.}
#'   }
#'
#' @return A list of class "forestsearch" containing:
#'   \describe{
#'     \item{sg.harm}{Character vector. Selected subgroup definition.}
#'     \item{grp.consistency}{List. Consistency evaluation results.}
#'     \item{find.grps}{List. Subgroup search results.}
#'     \item{confounders.candidate}{Character vector. Input confounders.}
#'     \item{confounders.evaluated}{Character vector. Evaluated confounders.}
#'     \item{df.est}{Data frame. Estimation dataset with subgroup flags.}
#'     \item{df.predict}{Data frame. Prediction dataset with subgroup flags.}
#'     \item{df.test}{Data frame. Test dataset with subgroup flags.}
#'     \item{minutes_all}{Numeric. Total computation time in minutes.}
#'     \item{grf_res}{List. GRF results (if use_grf = TRUE).}
#'     \item{args_call_all}{List. All function arguments for reproducibility.}
#'   }
