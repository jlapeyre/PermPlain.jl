# PermPlain

[![Build Status](https://travis-ci.org/jlapeyre/PermPlain.jl.svg?branch=master)](https://travis-ci.org/jlapeyre/PermPlain.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/PermPlain.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/permplain-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/jlapeyre/PermPlain.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/PermPlain.jl?branch=master)
[![codecov.io](http://codecov.io/github/jlapeyre/PermPlain.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/PermPlain.jl?branch=master)

This package implements functions for manipulating permutations.

See the docstrings `julia> ? PermPlain`.

The permutations are stored as

* Arrays of integers corresponding to one-line notation (representation)
* Arrays of arrays of integers corresponding to cycle notation (representation)
* A "sparse" indexable cycle notation

The methods do the following

* Generate permutations
* Convert between representations
* Compute properties of permutations
* Implement operations on and actions by permutations
* Print representations of permutations

The methods are meant to work easily with existing routines
for permutations in Base.

The cyclic decomposition is a represented by a array of arrays of Integers.
The array representation is an array of Integers.
If both input and output are permutations, then the input and output
representations are the same, unless otherwise noted.


## Some things to know

* The canonical form is: smallest element in a cycle is written first,
  cycles are sorted by increasing value of the first element.
