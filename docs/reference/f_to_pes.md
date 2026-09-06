# Convert Cohen's f to partial eta squared

Converts Cohen's f for an effect to its corresponding partial eta
squared: \$\$\eta_p^2 = \frac{f^2}{1 + f^2}.\$\$ The same formula
applies to between-subject, within-subject, and interaction effects when
f uses the error variance corresponding to that effect's partial eta
squared. This helper does not translate between G\*Power's
repeated-measures effect-size conventions.

## Usage

``` r
f_to_pes(f)
```

## Arguments

- f:

  A single finite, nonnegative numeric value representing Cohen's f, not
  an ANOVA F statistic or f-squared.

## Value

A single numeric partial eta squared. Zero maps to zero; power functions
retain their own restrictions on `target_pes`. Very large f values may
yield exactly one because of floating-point rounding.

## References

Cohen, J. (1988). Statistical power analysis for the behavioral sciences
(2nd ed.). Lawrence Erlbaum Associates.

## See also

[`compute_scale_factor()`](https://shaheedazaad.github.io/anovapowersim/reference/compute_scale_factor.md),
[`power_n_calc()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n_calc.md)

## Examples

``` r
f_to_pes(0.25) # approximately 0.05882
#> [1] 0.05882353
power_n_calc(
  between = c(group = 2),
  term = "group",
  target_pes = f_to_pes(0.25)
)
#> <anovapowersim_curve>
#>   term:          'group'
#>   target power:  0.900
#>   alpha:         0.05
#>   effect size:   pes = 0.0588
#>   n values:      13 per-cell sample sizes visited
#>   calculation:   calculated power only
#>   n needed for between-subjects cell: 87
#>   total N needed: 174
#> 
#>  n_per_cell total_n n_sims valid_sims failed_sims epsilon num_df den_df    ncp
#>           2       4     NA         NA          NA       1      1      2  0.125
#>           4       8     NA         NA          NA       1      1      6  0.375
#>           8      16     NA         NA          NA       1      1     14  0.875
#>          16      32     NA         NA          NA       1      1     30  1.875
#>          32      64     NA         NA          NA       1      1     62  3.875
#>          64     128     NA         NA          NA       1      1    126  7.875
#>          80     160     NA         NA          NA       1      1    158  9.875
#>          84     168     NA         NA          NA       1      1    166 10.375
#>          86     172     NA         NA          NA       1      1    170 10.625
#>          87     174     NA         NA          NA       1      1    172 10.750
#>          88     176     NA         NA          NA       1      1    174 10.875
#>          96     192     NA         NA          NA       1      1    190 11.875
#>         128     256     NA         NA          NA       1      1    254 15.875
#>  power_calc power_sim
#>       0.056      <NA>
#>       0.082      <NA>
#>       0.141      <NA>
#>       0.263      <NA>
#>       0.491      <NA>
#>       0.795      <NA>
#>       0.878      <NA>
#>       0.893      <NA>
#>       0.900      <NA>
#>       0.903      <NA>
#>       0.907      <NA>
#>       0.929      <NA>
#>       0.978      <NA>
```
