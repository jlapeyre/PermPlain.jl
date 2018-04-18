# PermPlain

[![Build Status](https://travis-ci.org/jlapeyre/PermPlain.jl.svg?branch=master)](https://travis-ci.org/jlapeyre/PermPlain.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/PermPlain.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/permplain-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/jlapeyre/PermPlain.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/PermPlain.jl?branch=master)
[![codecov.io](http://codecov.io/github/jlapeyre/PermPlain.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/PermPlain.jl?branch=master)

This package implements methods for manipulating permutations.
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
for permutations in Base. They are also wrapped by methods
on objects in the PermutationsA module.

The cyclic decomposition is a represented by a array of arrays of Integers.
The array representation is an array of Integers.
If both input and output are permutations, then the input and output
representations are the same, unless otherwise noted.

Usage:
```julia
p = randperm(10)          # permutation with Base.randperm
p = randperm(BigInt,10)   # added methods to select type of array elements.
p = randperm(Int32,10)    # 
c = permcycles(p)         # cyclic decompostion of p in canonical form.
cycstoperm(c)             # convert c to array form ( p == cycstoperm(c) is true )
cyclelengths(p)           # list of lengths of cycles in decomposition of p.
cyclelengths(c)           # list of lengths of cycles c.
cycletype(p)              # cycle type of p (as Accumulator).
cycletype(c)              # cycle type of c.
canoncycles(c)            # copy of c in canonical form.
permsgn(p)                # signature of p
permsgn(c)                # signature of c
permorder(p)              # order of permutation p
permorder(c)              # order of permutation c
permcommute(p,q)          # true if and only if p and q commute
permdistance(p,q)         # distance, ie number of points which have different image under p and q
permcompose(q,p)          # composition (multiplication) q * p
permcompose!(q,p)         # composition updating q
permpower(p,n)            # nth power of p
permpower(c,n)            # nth power of c
permkron(p,q)             # kronecker product of list permutations, induced by matrix kron product
cyc_pow_perm(c,n)         # nth power of c returned in array representation
m=permtomat(p,flag)       # permutation matrix, sparse if flag is true
mattoperm(m)              # array representation of permutation (abstract) matrix m
ipiv2perm(v,[n])          # convert pivot vector to permutation in list form
perm2piv(p)               # convert p to pivot vector
isperm(m)                 # true if matrix m is a permutation
isperm(c)                 # true if c is a permutation
Base.isperm(p)            # already present in Base module
permlistisequal(p,q)      # true if p and q represent the same permutation
                          # (the lists need not be of the same length)
ltpermlist(p,q)           # true if p < q lexicographically ( [1:n] is smallest permutation )
cycleprint(v)             # print vector v in form of a single cycle
permarrprint(p)           # print vector v in array (one-line) notation
copy(c)
copy(p)
```

A method ```isperm(A)``` that takes a matrix as input extends the ```Base.isperm()```.
```isperm``` returns true only if the matrix is a permutation matrix.

There are several othe methods that are not yet in this list.

## Some things to know

* Currently many method names are exported. But, I now prefer qualifying the names when calling the methods

* Cycles of length 1 are omitted from the disjoint-cycle and sparse representations

* The canonical form is: smallest element in a cycle is written first,
  cycles are sorted by increasing value of the first element.
