# Sampling a flat histogram along the cluster size #

## Files ##

 File         | Description
--------------|------------------------------
clusmc.c      | Monte Carlo, flat histogram,
clusmd.c      | molecular dynamics with hybrid MC, flat histogram



## clusmc Usage ##

./clusmc


## Output chist.dat ##

Column        | Meaning
--------------|------------------------------
1             | cluster size
2             | histogram (normalized)
3             | histogram (unnormalized)
4             | adaptive potential



## Notes ##

### Hybrid MC ###

The following notes apply to `clusmd.c`.

  * Higher acceptance ratios (> 50%) of hybrid MC is usually good, unlike regular MC.
  * Sometimes the system would be locked into a state in hybrid MC.  No solution when this happens.

#### Frequency of HMC ####
  * It appears to be most efficient to implement hybrid MC every single step.
  * Doing HMC less often decreases the acceptance ratio of HMC.

#### Velocity scrambling ####
  Scramble velocities after a (rejection of) hybrid MC step (see lj_vscramble() in ljcore.h).

  * One scramble is **needed** to achieve a flat histogram.
  * But more swaps do not help too much, it may be even counterproductive. 10 swaps are tried, counterproductive.
  * Swaps with 10% probability also hurts performance.

#### Time step of v-rescaling ####
  * Try to reduce the time step from 0.1 to 0.01.  No significant effect observed.

#### Initial lnf ####
  * Increasing lnf from 2e-5 to 1e-4 helps achieving a flat histogram faster.
  * Increasing lnf from 1e-4 to 1e-3 further helps, but the adaptive potential is messier.

#### Seed ####
  * Regularly changing the seed particle does not appear to have any effect.





## Ideas ##

* Superposition of ensembles of different parameters (e.g., densities)
  to the flat histogram.
  + ensembles should be well chosen
  + to avoid negative weight in combination

* Easy-hard coupling of the biased potential.
  X: hard order parameter, like RMSD
  Y: easy order parameter, like Go-energy
  Then the bias potential assumes the form of lambda(X) Y.

* Shape-based implicit Gibbs
  Draw a sphere and disable the interaction between
  particles inside and out of the sphere.
  Change the radius of the sphere.
  + currently, distribution focused on the two ends
  + missing a vol^n factor?
  + an elastic sphere, when changing the volume?
