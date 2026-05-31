# anovapowersim 1.1.0

* Added a tolerance argument to `power_n()` for more precise control over the adaptive search.

# anovapowersim 1.0.0

* First official release
* Fixed a bug where adaptive search for purely between-subjects designs would fail if the starting N was too small

# anovapowersim 0.2.0

* Added parallel processing for simulation runs in `power_curve()` and
  `power_n()`. Use `parallel = TRUE` to enable
  parallel simulations and `cores` to control the number of cores.
