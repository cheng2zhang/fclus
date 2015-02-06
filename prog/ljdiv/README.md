# Implicit Gibbs ensemble #

## Files ##

 File         | Description
--------------|------------------------------
igibbs.c      | implicit Gibbs ensemble
ljdiv.h       | division based energy (used by igibbs.c)


## Ideas ##

* Shape-based implicit Gibbs
  Draw a sphere and disable the interaction between
  particles inside and out of the sphere.
  Change the radius of the sphere.
  + currently, distribution focused on the two ends
  + missing a vol^n factor?
  + an elastic sphere, when changing the volume?
