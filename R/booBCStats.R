#' Differential connectivity bootstrap statistics (bias-corrected)
#'
#' Computes a two-sided p-value and confidence interval for a statistic \code{mi}
#' against its null distribution \code{pi}. Internally uses a bias-corrected
#' mapping on the null via \eqn{z_0} and finds the calibration \eqn{z_A} by
#' minimizing a CI width criterion (Brent's method).
#'
#' @param mi Numeric scalar. The observed measure/statistic.
#' @param pi Numeric vector. The null/sampling distribution (e.g., from
#'   permutation or bootstrap). \code{NA}s are removed.
#' @param conf Numeric in (0,1). Confidence level for the interval. Default 0.95.
#'
#' @return A one-row \code{data.frame} with columns:
#' \itemize{
#'   \item \code{test_dir}: "Up" if \code{mi > 0} else "Down"
#'   \item \code{z}: signed z-score aligned with \code{mi}
#'   \item \code{pval}: two-sided p-value
#'   \item \code{fdr}: Benjamini–Hochberg adjusted p-value (computed on the returned row)
#'   \item \code{ci_low}, \code{ci_high}: confidence interval bounds at level \code{conf}
#' }
#'
#' @details
#' Let \eqn{\hat{F}} be the empirical CDF of \code{pi}. Define
#' \deqn{ z_0 = \mathrm{sign}(-\mathbb{E}[pi - mi]) \cdot \left| \Phi^{-1}\{\hat{F}(mi)\} \right|, }
#' and find \eqn{z_A} that minimizes the CI width induced by the transformed
#' quantiles \eqn{\Phi(2 z_0 \pm z_A)}. The two-sided p-value is
#' \eqn{2 \Phi(-|z_A|)} and the (\code{conf})-level CI uses
#' \eqn{z_\alpha = \Phi^{-1}(\alpha)} with \eqn{\alpha = (1 - \mathrm{conf})/2}.
#'
#' @examples
#' set.seed(1)
#' mi <- 0.15
#' pi <- rnorm(1000, mean = 0, sd = 0.2)
#' bootBCStats(mi, pi)
#'
#' @importFrom stats pnorm qnorm quantile optimize p.adjust
#' @export
bootBCStats <- function(mi, pi, conf = 0.95) {
  # --- Validate inputs
  if (!is.numeric(mi) || length(mi) != 1L || is.na(mi)) {
    stop("`mi` must be a single, non-NA numeric value.")
  }
  if (!is.numeric(pi) || length(pi) < 10L) {
    stop("`pi` must be a numeric vector with length >= 10.")
  }
  pi <- pi[is.finite(pi)]
  if (length(pi) < 10L) {
    stop("After removing non-finite values, `pi` has < 10 observations.")
  }
  if (!is.numeric(conf) || length(conf) != 1L || conf <= 0 || conf >= 1) {
    stop("`conf` must be a single numeric in (0, 1).")
  }

  # --- Bias-correction term z0 (clip phat to avoid +/-Inf)
  phat_raw <- mean(pi <= mi)
  eps <- 1 / (length(pi) + 1)                      # continuity correction
  phat <- min(max(phat_raw, eps), 1 - eps)
  z0 <- abs(stats::qnorm(phat)) * sign(-mean(pi - mi))

  # --- Find zA by minimizing CI width induced by (2*z0 ± zA)
  obj <- function(zA) .getCI(za = zA, z0 = z0, pi = pi, return_interval = FALSE)
  zA  <- stats::optimize(f = obj, interval = c(-10, 10))$minimum

  # --- Two-sided p-value and signed z aligned with mi
  pval <- 2 * stats::pnorm(-abs(zA))
  z_signed <- abs(stats::qnorm(pval / 2)) * sign(mi)

  # --- Confidence interval at chosen level
  alpha <- (1 - conf) / 2
  z_alpha <- stats::qnorm(alpha)
  ci <- .getCI(za = z_alpha, z0 = z0, pi = pi, return_interval = TRUE)

  res <- data.frame(
    test_dir = ifelse(mi > 0, "Up", "Down"),
    z        = z_signed,
    pval     = pval,
    ci_low   = ci[1],
    ci_high  = ci[2]
  )
  res$fdr <- stats::p.adjust(res$pval, method = "BH")
  return(res)
}

#' @keywords internal
.getCI <- function(za, z0, pi, return_interval = FALSE) {
  q1 <- stats::pnorm(z0 * 2 + za)
  q2 <- stats::pnorm(z0 * 2 - za)
  if (return_interval) {
    stats::quantile(pi, sort(c(q1, q2)), names = FALSE, type = 7)
  } else {
    min(abs(stats::quantile(pi, c(q1, q2), names = FALSE, type = 7)))
  }
}
