# ExtremeVertexDesigns.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caseykneale.github.io/ExtremeVertexDesigns.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caseykneale.github.io/ExtremeVertexDesigns.jl/dev)
[![Build Status](https://github.com/caseykneale/ExtremeVertexDesigns.jl/workflows/CI/badge.svg)](https://github.com/caseykneale/ExtremeVertexDesigns.jl/actions)
[![Coverage](https://codecov.io/gh/caseykneale/ExtremeVertexDesigns.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/caseykneale/ExtremeVertexDesigns.jl)

## Status
Work in progress.

## Overview
Package for generating extreme value design candidate sets from inequality constraints. These types of experimental designs are useful for mixture experiments and generating calibrations.

## Example

```Julia
using ExtremeVertexDesigns
using JuMP

model = Model()

@variables(model, begin
    x
    y
end)

@constraints(model, begin
    0 <= x <= 1.0
    0 <= y <= 1.0
    x <= y
end)

evs = extreme_values( model; barycenters = true)

federov_D( 4, evs; k = 1e-6 )
```

## Credits
Thank you to Júlio Hoffimann, and Benoît Legat for helping me grok Polyhedra.jl and improving some code.