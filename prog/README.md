# Sampling a flat histogram along the cluster size #

## Files ##

 File         | Description
--------------|------------------------------
clusmc.c      | Monte Carlo, flat histogram,
clusmd.c      | molecular dynamics, flat histogram, using hybrid MC


## Notes ##

### Hybrid MC ###

  o Higher acceptance ratios (> 50%) of hybrid MC is usually good, unlike regular MC.
  o Sometimes the system would be locked into a state in hybrid MC.  No solution when this happens.

#### Frequency of HMC ####
  o It appears to be most efficient to implement hybrid MC every single step.
  o Doing HMC less often decreases the acceptance ratio of HMC.

#### Velocity scrambling ####
  Scramble velocities after a rejection of hybrid MC (see lj_vscramble() in ljcore.h).

  o One scramble is **needed** to achieve a flat histogram.
  o But more swaps do not help too much, it may be even counterproductive.

#### Time step of v-rescaling ####
  o Try to reduce the time step from 0.1 to 0.01.  No significant effect observed.

#### Initial lnf ####
  o increasing lnf from 2e-5 to 1e-4 helps achieving a flat histogram faster.
  o increasing lnf from 1e-4 to 1e-3 further helps, but the adaptive potential is messier.

