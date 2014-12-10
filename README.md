# PermPlain

[![Build Status](https://travis-ci.org/jlapeyre/PermPlain.jl.svg?branch=master)](https://travis-ci.org/jlapeyre/PermPlain.jl)

This package implements methods for manipulating permutations.
The methods operate on data types in the Base module, or in modules providing generic
data types. The permutations are stored as

* Arrays of integers corresponding to one-line notation (representation)
* Arrays of arrays of integers corresponding to cycle notation (representation)
* Matrices

The methods do the following

* Generate permutations
* Convert between representations
* Compute properties of permutations
* Implement operations on and actions by permutations
* Print representations of permutations

The methods are meant to work easily with existing routines
for permutations in Base. They are also wrapped by methods
on objects in the PermutationsA module.

This module is experimental and the interface should not
be considered stable.

Usage:
```julia
p = randperm(10)     # create permutation with Base.randperm
c = permcycles(p)    # compute cyclic decompostion of p
```

## Some things to know

* Cycles of length 1 are omitted from the disjoing cycle representation

* The canonical order is, smallest element in a cycle is written first,
  cycles are sorted by increasing value of the first element.
