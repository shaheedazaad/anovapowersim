# Move anovapowersim releases to its existing R-universe

## Summary

Keep the existing `shaheedazaad.r-universe.dev` page, replace its
automatic `main` tracking with tagged-release tracking, and publish
`1.2.0` as the final CRAN release carrying an attach-time migration
message.

## Package changes

- Add an `.onAttach()` startup message only when both conditions hold:

  - the installed version is exactly `1.2.0`;
  - the installed `Repository` metadata is `CRAN`.

- Display:

  > anovapowersim 1.2.0 is the final release on CRAN. Future versions
  > will be released through R-universe. Installation instructions:
  > <https://shaheedazaad.github.io/anovapowersim/#installation>

- Update the README/pkgdown Installation section to make R-universe the
  stable source:

  ``` r

  install.packages(
    "anovapowersim",
    repos = c(
      shaheedazaad = "https://shaheedazaad.r-universe.dev",
      CRAN = "https://cloud.r-project.org"
    )
  )
  ```

- Retain GitHub instructions for development snapshots and explain that
  CRAN remains frozen at `1.2.0`.

- Add the existing R-universe package page to `DESCRIPTION`, prepare
  version `1.2.0`, update `NEWS.md` and `cran-comments.md`, and rebuild
  the tracked pkgdown site.

No exported R API changes.

## Take control of the existing universe

- Confirm the R-universe GitHub app remains installed for the
  `shaheedazaad` account.

- Create `shaheedazaad/shaheedazaad.r-universe.dev` with:

  ``` json
  [
    {
      "package": "anovapowersim",
      "url": "https://github.com/shaheedazaad/anovapowersim",
      "branch": "*release"
    }
  ]
  ```

- Let this custom registry replace the current automatically generated
  registry, which tracks `HEAD`.

- Tag the final commit as `v1.2.0` and create a GitHub Release. Future
  stable versions will use the same tag-and-release process.

## Verification and rollout

- Test that the message appears only for CRAN-installed `1.2.0`; confirm
  silence for development, R-universe, GitHub, local, and later
  versions.
- Run the test suite and `R CMD check`, rebuild pkgdown, and verify the
  Installation link.
- Confirm R-universe reports `RemoteRef` as the latest release rather
  than `HEAD`.
- Install `1.2.0` from R-universe in a temporary library and confirm
  quiet attachment.
- Submit `1.2.0` to CRAN, then install its published build and confirm
  the migration message appears.
