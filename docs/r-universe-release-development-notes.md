# R-universe release and development channels

Checked 2026-09-09 against official documentation.

The existing `migration_plan.md` already chooses a workable separation: stable
releases on `shaheedazaad.r-universe.dev`, development snapshots installed from
GitHub.

R-universe supports one version of each package per universe. The registry's
`branch` field accepts a Git branch or tag; `"*release"` follows the latest
GitHub release. Omitting `branch` tracks the default branch. Consequently, the
plan's `"branch": "*release"` prevents ordinary development commits on `main`
from becoming the published stable package. Create a GitHub Release for each
stable version; a Git tag alone is not the release workflow described by this
setting. An explicit tag such as `"v1.2.0"` is an alternative when manual
promotion is preferred. [R-universe setup documentation](https://docs.r-universe.dev/publish/set-up.html#tracking-custom-branches-or-releases)

Recommended development installation:

```r
remotes::install_github("shaheedazaad/anovapowersim", ref = "main")
```

The `ref` argument selects a branch, tag, or commit; its default is the
repository's default branch. This builds from GitHub source, so development
users do not receive R-universe's prebuilt package binaries.
[remotes documentation](https://remotes.r-lib.org/reference/install_github.html)

If both channels need R-universe binaries, they need separate universes. The
official documentation ties a universe to a GitHub account/organization and
describes separate organizations for development and production. Both registries
can point at the same package repository, with different `branch` values; the
source repository does not need to move. This is additional administration and
is unnecessary for the current plan.
[R-universe setup documentation](https://docs.r-universe.dev/publish/set-up.html#can-a-user-have-multiple-universes-or-multiple-versions-of-a-package)

When offering two repository channels, keep development out of the stable
installation command: R's default duplicate-package filtering selects the
highest version across repositories, using repository order only for equal
versions. [R available.packages documentation](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/available.packages.html)

An optional future alternative is R-multiverse Community for published releases
and the personal universe for development. Community requires registration and
then deploys source-code releases automatically. Its separate Production
repository adds periodic snapshots and quality requirements, so these two
R-multiverse services should not be conflated.
[Community](https://r-multiverse.org/community.html),
[Production](https://r-multiverse.org/production.html)
