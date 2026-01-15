# Differential connectivity bootstrap statistics (bias-corrected)

Computes a two-sided p-value and confidence interval for a statistic
`mi` against its null distribution `pi`. Internally uses a
bias-corrected mapping on the null via \\z_0\\ and finds the calibration
\\z_A\\ by minimizing a CI width criterion (Brent's method).

## Usage

``` r
.bootBCStats(mi, pi, CI = 0.95)
```

## Arguments

- mi:

  Numeric scalar. The observed measure/statistic.

- pi:

  Numeric vector. The null/sampling distribution (e.g., from permutation
  or bootstrap). `NA`s are removed.

- CI:

  Numeric in (0,1). Confidence level for the interval. Default 0.95.

## Value

A one-row `data.frame` with columns:

- `test_dir`: "Up" if `mi > 0` else "Down"

- `z`: signed z-score aligned with `mi`

- `pval`: two-sided p-value

- `fdr`: Benjamini–Hochberg adjusted p-value (computed on the returned
  row)

- `ci_low`, `ci_high`: confidence interval bounds at level `conf`

## Details

Let \\\hat{F}\\ be the empirical CDF of `pi`. Define \$\$ z_0 =
\mathrm{sign}(-\mathbb{E}\[pi - mi\]) \cdot \left\|
\Phi^{-1}\\\hat{F}(mi)\\ \right\|, \$\$ and find \\z_A\\ that minimizes
the CI width induced by the transformed quantiles \\\Phi(2 z_0 \pm
z_A)\\. The two-sided p-value is \\2 \Phi(-\|z_A\|)\\ and the
(`conf`)-level CI uses \\z\_\alpha = \Phi^{-1}(\alpha)\\ with \\\alpha =
(1 - \mathrm{conf})/2\\.
