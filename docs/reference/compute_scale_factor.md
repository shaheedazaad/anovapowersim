# Compute the mean-deviation scaling factor from a change in partial eta squared

Given an existing partial eta squared for a term and a target partial
eta squared, returns the multiplier `k` that must be applied to that
term's additive contribution to the cell means in order to obtain the
target effect size under the same residual structure.

## Usage

``` r
compute_scale_factor(old_pes, new_pes)
```

## Arguments

- old_pes:

  Numeric scalar in (0, 1), or a numeric-looking character scalar such
  as `".310"`. The current partial eta squared for the term of interest.

- new_pes:

  Numeric scalar in (0, 1), or a numeric-looking character scalar such
  as `".200"`. The target partial eta squared.

## Value

A single positive numeric value `k`. `k > 1` amplifies the effect,
`k < 1` shrinks it, and `k == 1` leaves it unchanged.

## Details

The derivation is straightforward: partial eta squared can be written as
`pes = SS_effect / (SS_effect + C)`, where `C` is the part of the
denominator held fixed by this package's rescaling. Thus
`pes / (1 - pes)` scales as the target effect's sum of squares. Scaling
the term's deviations by `k` scales the target effect's sum of squares
by `k^2`, so the required multiplier is \$\$k =
\sqrt{\frac{p\_{\mathrm{new}} / (1 - p\_{\mathrm{new}})}
{p\_{\mathrm{old}} / (1 - p\_{\mathrm{old}})}}.\$\$

## See also

[`design_term_means()`](https://shaheedazaad.github.io/anovapowersim/reference/design_term_means.md),
[`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)

## Examples

``` r
compute_scale_factor(0.10, 0.05)   # shrink
compute_scale_factor(0.05, 0.10)   # amplify
```
