# convert_length_to_mass

A function to convert species-specific lengths to mass based on a
user-provided length-mass equation form and variable values

## Usage

``` r
convert_length_to_mass(
  taxaSampleList = NULL,
  taxaInfo = NULL,
  lengthValue = NULL,
  massValue = NULL,
  reduce = TRUE,
  ...
)
```

## Arguments

- taxaSampleList:

  data.frame in long format for a single taxa. The data.frame should
  contain a species identifier column \`taxonID\` and a column of length
  bin categories lengthValue.

- taxaInfo:

  a data.frame of the information to convert length to mass for all
  taxa. The taxa specified in \`taxaSampleList\` will be subset from
  here. This data.frame must contain a \`taxonID\` column, the
  length-to-mass equation formula, \`massForm\`, which must contain
  \`lengthClass\` as a variable (e.g., \`afdm_mg~a\*lengthClass^b\`).
  Additional columns are necessary based on the length-mass formula. All
  other non-\`lengthClass\` variables on the right hand side (RHS) must
  have unique columns named the variable name. For example, the above
  formula structure, \`afdm_mg~a\*lengthClass^b\`, the RHS is
  \`a\*lengthClass^b\`. \`lengthClass\` is a required column, but
  optionally necessary columns are \`a\` and \`b\` for the other
  variables. The formula will be parsed and species-specific \`a\` and
  \`b\` coefficients will be inserted for conversion.

- lengthValue:

  a character string identifying the column name for the length class.
  The data should be numeric or coercible

- massValue:

  a character string to name the mass variable. This should not contain
  any punctuation (e.g., .,\_,-, etc.).

- reduce:

  logical. If TRUE (default) the mass column will be added to
  \`taxaSampleList\`. The name of the mass column will be parsed from
  the left hand side (LHS) of the \`massForm\` provided in \`taxaInfo\`.
  For example, in \`mass~a\*lengthValue^b\` the mass column will be
  named based on the massValue column. This data.frame is returned if
  reduce == TRUE

- ...:

  additional arguments passed to function

## Value

taxaSampleList with the mass column added.
