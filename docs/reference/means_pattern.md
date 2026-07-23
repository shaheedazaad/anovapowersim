# Define a sparse relative cell-mean pattern

Creates a sparse mean-shape specification for the balanced simulation
functions. End each cell definition with `value`. Unlisted cells have
raw value zero, and factors in the tested term that are omitted from
every row are broadcast when the pattern is resolved against a design.

## Usage

``` r
means_pattern(...)
```

## Arguments

- ...:

  Repeated named sparse-cell definitions. Each definition must use the
  same factor names in the same order and end in a finite numeric scalar
  named `value`. Factor levels may be supplied as one-based integer
  indices (for example, `time = 3`) or as the generated balanced-design
  names (for example, `time = "time3"`). The two forms are equivalent
  after the pattern is resolved against a design.

## Value

An object of class `anovapowersim_means_pattern`, retaining the sparse
definitions until a balanced simulation function resolves them against
its design and tested term.

## Details

Pattern values describe relative shape, not effect magnitude. The
selected power function projects the raw values onto the requested ANOVA
term, normalizes that component, and rescales it uniformly to reach
`target_pes`. Multiplying all values by one positive constant, adding an
intercept or a lower-order component, or reversing every sign therefore
leaves the same target-term direction (up to sign). Under nonsphericity,
different directions within a multi-df term can nevertheless produce
different simulated power.

This differs from
[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md),
where each `m` is a literal population mean whose magnitude directly
determines the simulated effect.

## Default direction

When no pattern is supplied, balanced simulations use centered scores in
generated level order, `i - (L + 1) / 2` for levels `i = 1, ..., L`,
normalized to unit length. Interactions use the Kronecker product of
their factors' normalized score vectors, followed by one final
normalization after broadcasting. This is an ordered, reproducible
convention rather than a neutral scientific assumption; an explicit
pattern is recommended whenever the expected shape is known.

## See also

[`cell_design()`](https://shaheedazaad.github.io/anovapowersim/reference/cell_design.md)
for unbalanced designs with literal cell means.

## Examples

``` r
trend <- means_pattern(
  time = 1, value = 0,
  time = 2, value = 0.3,
  time = 3, value = 0.5,
  time = 4, value = 0.6
)

interaction_shape <- means_pattern(
  group = "group1", time = "time3", value = 1,
  group = "group2", time = "time3", value = -1
)
```
