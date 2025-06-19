# cliffs-src

Source code for a numerical model (Cliffs) to simulate tsunami propagation and run-up on land.

### Cliffs Features

- Shallow-Water approximation with an option to emulate physical dispersion
- Use of Cartesian or spherical (lon/lat) coordinates
- 1D and 2D configurations, including 1-D channels with varying width - domain type is detected automatically
- Structured colocated grid with (optionally) varying spacing
- Runup on land
- Variety of starting/forcing options: bottom deformation, surface deformation, boundary input, or combination of those
- Output options: various
- Grid nesting with one-way coupling
- Parallelized with OpenMP
- NetCDF format of input/output data

### Modeling set-up examples

- [*Runup onto a sloping beach*](Examples/slpbeach.zip)
- [*Hokkaido 1993 (Okushiri) Tsunami*](Examples/OkushiriTsunami.zip)

(to run the examples, see [*Cliffs User Manual*](http://arxiv.org/abs/1410.0753) for instructions) 

### Documentation and References

- E. Tolkova. [*Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting*](http://link.springer.com/10.1007/s00024-014-0825-8) Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314
- [*Cliffs User Manual*](http://arxiv.org/abs/1410.0753)
- [*Cliffs Benchmarking Results*](https://arxiv.org/abs/1601.06486)
- [*Comparative simulations of the 2011 Tohoku tsunami with MOST and Cliffs*](https://arxiv.org/abs/1401.2700)