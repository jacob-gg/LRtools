#' \code{spline_refit()}
#'
#' @description \code{spline_refit()} takes an \code{lm} or \code{glm} model
#' object (generally containing linear terms) and refits it with nonlinear
#' spline terms. The function returns the refit model (invisibly) and prints
#' console output with AIC and BIC values from the original and new models
#' for comparison.
#'
#' @details If the user leaves the \code{to_spline} argument at its default
#' value, \code{"all"}, the function tries to identify main effect terms that
#' are specified as linear in the model and refit them as splines.
#' If \code{to_spline} is \code{"all"}, then the number given as the \code{df}
#' argument applies to all new spline terms. Alternatively, the user can
#' specify a set of linear main effects to refit as splines via the
#' \code{to_spline} argument. In that case, the user can set \code{df} to be a
#' single value (which will serve as the df of all new spline terms), or the
#' user can pass a vector of df values to \code{df}, with each corresponding to
#' a desired spline term.
#'
#' Currently, the function will only refit with splines terms that can be
#' identified as straightforward main effects of continuous variables:
#' It won't refit interaction terms (if created via \code{y ~ x*z} or
#' \code{y ~ x:z}), terms for binary predictors or factor variables,
#' terms created with \code{I()}, or any other terms containing \code{:}
#' (signaling an interaction) or parentheses (usually signaling an on-the-fly
#' variable computation or rescale). The primary limitation of this is that if
#' a user wants to fit a spline term to a predictor that has been put on the
#' log scale, they must create the logged version of the variable in the data
#' set underlying the model and use that variable directly in the model.
#'
#' @param mod A model object that inherits from \code{lm} or \code{glm}.
#' @param to_spline The terms to refit as spline terms (\code{"all"} or a
#' character vector).
#' @param df Degrees of freedom for the spline term(s) (either a single value,
#' in which case all spline terms will be fit with that df, or a numeric vector
#' with the same length as \code{to_spline}).
#' @param spline_type Use cubic or natural (restricted cubic) splines?
#' Must be either \code{"bs"} or \code{"ns"}.
#'
#' @return The updated model (invisibly), plus output in the console.
#'
#' @examples
#' m <- lm(mpg ~ hp, mtcars)
#' spline_refit(mod = m, to_spline = 'all', df = 4, spline_type = 'ns')
#'
#' # Function identifies the binary predictor (am) and the
#' # interaction term and only refits drat with splines:
#' m <- glm(vs ~ drat*am, mtcars, family = 'binomial')
#' spline_refit(m, 'all', 4, 'ns')
#'
#' # To refit the logged wt term in the model below with splines,
#' # one would need to compute log(wt) in the data set and use the
#' # corresponding variable in the model
#' m <- lm(log(mpg) ~ hp + wt + log(wt), mtcars)
#' spline_refit(m, c('hp', 'wt'), 4, 'bs')
#'
#' # The updated model is returned invisibly:
#' m <- glm(vs ~ drat*am, mtcars, family = 'binomial')
#' updated_mod <- spline_refit(m, 'all', 4, 'ns')
#' ggeffects::ggpredict(updated_mod, terms = c('drat', 'am')) |>
#'   plot()
#'
#' @export
spline_refit <- function(mod, to_spline = 'all', df = 5, spline_type = 'ns') {
  stopifnot('mod must be lm or glm class' =
              inherits(mod, 'lm') || inherits(mod, 'glm'),
            'to_spline must be "all" or a vector of (main effect) terms present in mod' =
              (length(to_spline) == 1 && to_spline == 'all') ||
              (length(to_spline) >= 1 & all(to_spline %in% attr(mod$terms, 'term.labels'))),
            'df must be a single number or a list of values equal in length to to_spline' =
              (length(df) == 1 || (length(df) == length(to_spline))) && is.numeric(df),
            'spline_type must be "ns" or "bs"' =
              length(spline_type) > 0 && spline_type %in% c('ns', 'bs'))

  # Get terms from input model and separate apparent main effects of
  # non-binary, non-factor variables from all other terms (interactions,
  # binary predictors, factors, I() terms, etc, e.g., x:z, I(x/z)
  initial_terms <- attr(mod$terms, 'term.labels')
  apparent_continuous_mains <- initial_terms[grepl('\\:|\\)|\\(', initial_terms) == F]
  binary_or_factor <- sapply(mod$model[apparent_continuous_mains],
                             \(x) ifelse(length(unique(x)) == 2 || is.factor(x),
                                         TRUE, FALSE))
  apparent_continuous_mains <- apparent_continuous_mains[!binary_or_factor]
  all_other_terms <- initial_terms[!(initial_terms %in% apparent_continuous_mains)]

  # Generate term set for updated model, now including splines as requested
  if (length(to_spline) == 1 && to_spline == 'all') {
    new_terms <- c(paste0('splines::', ifelse(spline_type == 'ns', 'ns(', 'bs('),
                          apparent_continuous_mains, ', df = ', df, ')'),
                   all_other_terms)
  } else {
    new_terms <- c(paste0('splines::', ifelse(spline_type == 'ns', 'ns(', 'bs('),
                          apparent_continuous_mains[apparent_continuous_mains %in% to_spline], ', df = ', df, ')'),
                   apparent_continuous_mains[apparent_continuous_mains %in% to_spline == FALSE],
                   all_other_terms)
  }

  # Fit updated model
  mod_out <- stats::update(mod, stats::reformulate(new_terms, stats::formula(mod)[[2]]))

  # Prep and print comparative model info
  aic <- stats::AIC(mod, mod_out)
  bic <- stats::BIC(mod, mod_out)
  # Get max nchars of all 2-digit rounded AIC values, and pick highest plus 1
  # (or 4, if 4 is greater) as number of characters to print per AIC "cell"
  prnt_lngth_aic <- max(4, max(nchar(sprintf('%.2f', round(aic$AIC, 2)))) + 1)
  # Pick max nchar of all integer model df values, and pick highest plus 1
  # (or 3, if 3 is greater) as number of characters to print per df "cell"
  prnt_lngth_df <- max(3, max(nchar(aic$df)) + 1)
  # Format and print
  mod_formulas_ln1 <- paste0('\n\rOld model formula: ',
                             deparse(stats::formula(mod), width.cutoff = 500))
  mod_formulas_ln2 <- paste0('\n\rNew model formula: ',
                             deparse(stats::formula(mod_out), width.cutoff = 500))
  glm_info <- ifelse(inherits(mod_out, 'glm'),
                     paste0('\n\n\rFamily: ', mod_out$family$family,
                            '\n\rLink: ', mod_out$family$link,
                            '\n'),
                     '\n')
  info_criteria_ln1 <- sprintf('\n%-11s%-*sAIC%-*sBIC',
                               '', prnt_lngth_df, 'df', prnt_lngth_aic - 3, '')
  info_criteria_ln2 <- sprintf('\nOld model: %-*d%-*.2f%-.2f',
                               prnt_lngth_df, aic$df[1],
                               prnt_lngth_aic, aic$AIC[1], bic$BIC[1])
  info_criteria_ln3 <- sprintf('\nNew model: %-*d%-*.2f%-.2f',
                               prnt_lngth_df, aic$df[2],
                               prnt_lngth_aic, aic$AIC[2], bic$BIC[2])
  cat(mod_formulas_ln1,
      mod_formulas_ln2,
      glm_info,
      info_criteria_ln1,
      info_criteria_ln2,
      info_criteria_ln3,
      '\n')

  invisible(mod_out)
}
