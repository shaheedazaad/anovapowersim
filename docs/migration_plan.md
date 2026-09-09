# Track releases on R-universe while retaining CRAN

## Current decision

The proposed migration away from CRAN is cancelled. CRAN remains the
primary installation source, with continuing CRAN releases. There is no
final-CRAN version or attach-time migration message.

R-universe is an additional source for published GitHub Releases, not
development snapshots. Development versions remain available from
GitHub.

## R-universe configuration

The registry at `shaheedazaad/shaheedazaad.r-universe.dev` contains:

``` json
[
  {
    "package": "anovapowersim",
    "url": "https://github.com/shaheedazaad/anovapowersim",
    "branch": "*release"
  }
]
```

Publish a GitHub Release for each version intended for R-universe.
Ordinary commits to `main` do not change its published version. CRAN
submission remains an independent release step; publishing a GitHub
Release does not submit to CRAN.

## Installation

``` r

install.packages("anovapowersim")
```

For development snapshots:

``` r

remotes::install_github("shaheedazaad/anovapowersim")
```
