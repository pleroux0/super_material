# TODO List

This list is here to keep track of things that can still be done, but not things that have to be done.

### Mattis-Bardeen superconductor conductivity

The Mattis-Bardeen superconductor conductivity implementation is fully functional. Only numerical stability, performance, and validation improvements are left.

* Analytically perform all the integral interval transforms
* Perform common subexpression elimination on integrand evaluations
* Find out where the sign inversion in the Mattis-Bardeen expression comes from

### Zimmermann superconductor conductivity

The Zimmermann superconductor conductivity implementation is fully functional.  Only numerical stability and performance improvements are left.

* Analytically perform all the integral interval transforms
* Find improved integral interval transforms
* Perform common subexpression elimination on integrand evaluations

### Gap Energy

* Add `BCSWeakGapEnergy` for faster weak-coupling calculations

### Documentation

The public API is fully documented, but there are improvements that can be made.

* Move more of the documentation to be inline
* Changes necessary to use apidoc to generate API documentation
* Link classes in type hints
* Do not generate fully qualified names in parameter types
