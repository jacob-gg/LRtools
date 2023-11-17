#' \code{litl()}
#'
#' @description \code{litl()} = \code{l}inear \code{i}n \code{t}he \code{l}ogit.
#' \code{litl()} generates plots for graphically assessing the reasonableness
#' of the linearity assumption for one- and two-variable logistic regression
#' models under consideration.
#'
#' @details \code{litl()} assesses linearity in the logit between one or two
#' predictor variables, \code{x} and \code{z}, and a binary outcome, \code{y}.
#' \code{x} and \code{y} are required: \code{litl()} divides \code{x} into
#' \code{ntiles_x} quantile buckets (10 by default), calculates the log-odds
#' of \code{y} in each bucket, and generates a corresponding plot, as described
#' \href{https://library.virginia.edu/data/articles/graphical-linearity-assessment-one-and-two-predictor-logistic-regressions}{here}
#' and in Harrell (2015, \emph{Regression Modeling Strategies}).
#' If \code{z} is given, observations are stratified before quantile groups are
#' formed based on \code{x}:
#'
#' \itemize{
#'  \item{If \code{z} is numeric and not binary, and if a numeric value of
#'  \code{ntiles_z} is given, observations are stratified into \code{ntiles_z}
#'  quantile groups based on \code{z}, and observations within groups are then
#'  bucketed based on \code{x}.}
#'  \item{If \code{z} is not numeric and/or if \code{z} is binary, observations
#'  are stratified based on unique values of \code{z} before they are bucketed
#'  based on \code{x}.}}
#'
#'  By default, a least-squares line(s) of best fit is added to the resulting
#'  plot: If no \code{z} is provided, one line is added across buckets of
#'  \code{x}; if \code{z} is provided, a line is added for each strata of
#'  \code{z}. A LOESS line can be added instead by setting \code{plot_line =
#'  'loess'}, and the line(s) can be omitted by setting \code{plot_line = NULL}.
#'
#' @param x A double, integer, or logical predictor variable.
#' @param y A logical (\code{TRUE/FALSE/NA}) or binary numeric (\code{1/0/NA})
#' outcome variable.
#' @param z An optional second predictor variable, handled variously depending
#' on its type (see \bold{Details}).
#' @param ntiles_x A double or integer value indicating how many buckets to
#' divide \code{x} into (10 by default).
#' @param ntiles_z An optional double or integer value indicating how many
#' buckets to divide \code{z} into (if \code{z}) is provided.
#' @param plot_line The type of line to be fit in the linearity-assessment plot:
#' \code{'lm'} (default), \code{'loess'}, or \code{NULL} (no line).
#'
#' @return Invisibly, a two-element list: The first element, \code{[['data']]},
#' contains the data displayed in the plot; the second, \code{[['plot']]}, is
#' the plot object, with class \code{ggplot}.
#'
#' \code{x}, \code{y}, and \code{z} (if given) must be the same length.
#' \code{NA} values are allowed, although they are removed before the
#' calculative work proceeds.
#'
#' If any of the buckets of \code{x} (optionally stratified by \code{z}) have
#' < 20 cases, a warning is thrown. If P(Y=1) is 1 or 0 in any bucket(s),
#' rendering log(P/(1-P)) undefined, a note indicating so (and indicating where)
#' appears on the plot.
#'
#' @examples
#' set.seed(11)
#' xx <- rnorm(n = 2500, mean = 1, sd = 1)
#' zz <- rbinom(2500, 1, .5)
#' logit_y <- .5*xx + 1*zz + .5*xx*zz
#' prob_y <- exp(logit_y) / (1 + exp(logit_y))
#' yy <- rbinom(length(prob_y), 1, prob_y)
#' yy <- ifelse(yy == 1, TRUE, FALSE)
#' litl(x = xx, y = yy, z = zz)
#'
#' @importFrom rlang .data
#'
#' @export
litl <- function(x, y, z = NULL, ntiles_x = 10, ntiles_z = NULL, plot_line = 'lm') {
  # Validate input
  if (typeof(x) %in% c('double', 'integer', 'logical') == FALSE) {
    stop('x must be double, integer, or logical', call. = FALSE)
  }
  if ((typeof(y) == 'logical' || (typeof(y) %in% c('double', 'integer') && all(unique(y) %in% c(0, 1, NA)))) == FALSE) {
    stop('y must be logical (T/F/NA) or binary numeric (1/0/NA)', call. = FALSE)
  }
  if (length(x) != length(y)) {
    stop('x and y must be the same length', call. = FALSE)
  }
  if (typeof(ntiles_x) %in% c('double', 'integer') == FALSE) {
    stop('ntiles_x must be double or integer', call. = FALSE)
  }
  if (plot_line %in% c('lm', 'loess') == FALSE && is.null(plot_line) == FALSE) {
    stop('plot_line must be "lm", "loess", or NULL', call. = FALSE)
  }

  # Structure data for later splitting and recombining
  if (is.null(z)) {
    d <- data.frame(x, y, z = 1)
  } else {
    if (length(z) != length(x)) {
      stop ('z, x, and y must be the same length', call. = FALSE)
    }
    # If z is numeric and not binary, chop it up into ntiles_z buckets (although
    # if ntiles_z is not valid, ignore, warn, and use terciles instead); if z is
    # not numeric and/or is binary, let split() later divide it by unique values
    if ((typeof(z) %in% c('double', 'integer')) && (all(unique(z) %in% c(0, 1, NA)) == FALSE)) {
      if (is.null(ntiles_z) == FALSE && (typeof(ntiles_z) %in% c('double', 'integer'))) {
        z <- cut(z,
                 breaks = stats::quantile(z, probs = seq(0, 1, length = ntiles_z + 1), na.rm = TRUE),
                 include.lowest = T)
      } else {
        warning('ntiles_z is NULL or non-numeric: dividing z by terciles', call. = FALSE)
        z <- cut(z,
                 breaks = stats::quantile(z, probs = seq(0, 1, length = 3 + 1), na.rm = TRUE),
                 include.lowest = T)
      }
    }
    d <- data.frame(x, y, z)
  }

  # Check for NAs and report removals
  if (anyNA(d)) {
    d <- d[stats::complete.cases(d), ]
    warning(paste0('input variables contained NAs: dropping ', length(x) - nrow(d), ' incomplete cases'), call. = FALSE)
  }

  # Stratify data by z (if binary/categorical) or by quantiles of z (if z is numeric and valid ntiles_z provided)
  strat_dat <- split(d, d$z)

  # Determine quantiles of x for each z-stratified chunk of cases, then recombine strata
  for (i in 1:length(strat_dat)) {
    strat_dat[[i]]$ntiles <- cut(strat_dat[[i]]$x,
                                 breaks = stats::quantile(strat_dat[[i]]$x,
                                                          probs = seq(0, 1, length = ntiles_x + 1)),
                                 include.lowest = T)
  }
  d <- do.call('rbind', strat_dat)

  # Check N within quantiles and report if counts are < 20
  bucket_ns <- stats::aggregate(d$y, by = list(d$z, d$ntiles), function(n) length(n))$x
  if (any(sum(bucket_ns < 20))) {
    warning(paste0(sum(bucket_ns < 20),
                   ' ntiles(s) contain(s) < 20 observations\n',
                   'consider reducing ntiles_x and/or ntiles_z\n',
                   'log-odds estimates may be very noisy'),
            call. = FALSE)
  }

  # Structure data for plotting
  to_plt <- stats::aggregate(d$y, by = list(d$z, d$ntile), mean)
  colnames(to_plt) <- c('z', 'ntile_x', 'prob')
  to_plt$log_odds <- log(to_plt$prob / (1 - to_plt$prob))
  to_plt$ntile_x <- rep(1:ntiles_x, times = length(strat_dat))
  to_plt$z <- factor(to_plt$z)

  # Conditionally initialize one or another plot depending on whether z-stratification occurred
  if (length(unique(z)) > 1) {
    p <- ggplot2::ggplot(to_plt, ggplot2::aes(.data$ntile_x, .data$log_odds, color = z)) +
      ggplot2::labs(x = 'n-tile of x', y = 'Log-odds of y = 1', color = 'z',
                    title = 'Estimated log-odds of y = 1 across n-tiles of x, stratified by z')
  } else {
    p <- ggplot2::ggplot(to_plt, ggplot2::aes(.data$ntile_x, .data$log_odds)) +
      ggplot2::labs(x = 'n-tile of x', y = 'Log-odds of y = 1',
                    title = 'Estimated log-odds of y = 1 across n-tiles of x')
  }

  # Construct rest of plot
  p <- p +
    ggplot2::geom_point() +
    ggplot2::geom_line(alpha = .25) +
    { if (is.null(plot_line) == FALSE) ggplot2::geom_smooth(method = plot_line, se = FALSE, lwd = .5) } +
    ggplot2::scale_x_continuous(breaks = 1:ntiles_x) +
    ggplot2::geom_text(ggplot2::aes(.data$ntile_x, .data$log_odds,
                                    label = ifelse(is.infinite(.data$log_odds),
                                                   paste0('Undefined\nP(Y=1) = ',
                                                          ifelse(.data$log_odds > 0, '1', '0')),
                                                   ''),
                                    vjust = ifelse(is.infinite(.data$log_odds),
                                                   ifelse(.data$log_odds > 0, 2, -2),
                                                   0)),
                       size = 2.5, color = 'black') +
    ggplot2::theme(legend.position = 'bottom')

  suppressWarnings(suppressMessages(print(p)))
  invisible(list(data = to_plt, plot = p))
}
