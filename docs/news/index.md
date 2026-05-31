# Changelog

## anovapowersim 1.1.0

CRAN release: 2026-05-31

- Added a tolerance argument to
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md)
  for more precise control over the adaptive search.

## anovapowersim 1.0.0

CRAN release: 2026-05-28

- First official release
- Fixed a bug where adaptive search for purely between-subjects designs
  would fail if the starting N was too small

## anovapowersim 0.2.0

- Added parallel processing for simulation runs in
  [`power_curve()`](https://shaheedazaad.github.io/anovapowersim/reference/power_curve.md)
  and
  [`power_n()`](https://shaheedazaad.github.io/anovapowersim/reference/power_n.md).
  Use `parallel = TRUE` to enable parallel simulations and `cores` to
  control the number of cores.
