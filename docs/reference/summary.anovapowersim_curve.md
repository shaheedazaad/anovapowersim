# Summarise an anovapowersim power curve

Returns the full `$results` tibble along with a small header containing
the target, effective effect size, and estimated `n_needed`.

## Usage

``` r
# S3 method for class 'anovapowersim_curve'
summary(object, ...)
```

## Arguments

- object:

  An `anovapowersim_curve` object.

- ...:

  Unused.

## Value

A list with elements `header` (named character) and `curve` (tibble),
invisibly; printed to console as well.
