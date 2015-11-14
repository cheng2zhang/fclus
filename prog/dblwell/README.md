# Overview


This directory contains a program `dblwell`
which is intended to show the application of
hybrid Monte Carlo (HMC), a.k.a. Hamiltonian Monte Carlo,
in sampling a distribution with discontinuity.


# Theory

The potential in question is
$$
U(x) = (1/2) (x - y(x))^2,
$$
where $y(x)$ is $1.0$ if $x > xc$ or $-1.0$ otherwise
(the default $xc$ is $-0.5$).
Thus, the potential is discontinuous at $x = xc$.

If we ignore the $x$ dependence in $y$,
the approximate force is given by
$$
\tilde f = -\partial U/\partial x.
$$
Using this approximate force to guide the MD
would yield a distribution that is continuous at $x = xc$.
We can use HMC to fix the situation.
The idea is that we can regularly reject a piece of trajectory
that is about to cross the discontinuity,
and thence achieve the intended distribution.

Mathematically, after one or a few MD steps, we shall compute
a quantity $Q$.  Then we decide whether to accept or
to reject this short piece of trajectory, according to the probability
$$
A = \min[1, \exp(-\beta Q)].
$$
When a rejection occurs, we have to revert the velocity as well.



## Computing the quantity $Q$

There are three ways to compute the quantities.

1. The first is a formula that compute the total energy change $Q = \Delta E$.
If the exact force were used, the total energy should be a conserved.
The excess is due to the error of the propagation, and thus, should be deducted.
We refer to this formula the total energy formula.
Note that in using this formula, we only include the total energy change
caused by the MD integrator, and the change caused by the thermostat is excluded.

2. The second formula is derived from a transformation the above quantity
by Hamilton's equations
$$
dE/dt = (\partial E/\partial p) (dp/dt) + (\partial E/\partial x) (dx/dt)
      + (\partial E/\partial y) (dy/dt),
$$
where $p = m v$ denote the momentum.
By using the approximate force the first two terms cancel, and
$$
\Delta E = (\partial E/\partial y) \Delta y = (\partial U/\partial y) \Delta y.
$$
In practice, the energy change in each time step
can be computed as the average of the value at the beginning
and that at the end of the step.
$$
\Delta E
\approx (1/2) [ (\partial U/\partial y)_{x1}
              + (\partial U/\partial y)_{x2} ] \Delta y
\approx (1/2) [ U(x1, y2) - U(x1, y1) + U(x2, y2) - U(x2, y1) ].
$$
where $x1$ and $x2$ denote the coordinates at the beginning
and the end of an MD step, respectively.
The advantage of this method is that it depends only on
the coordinates.
This is referred to as the four-potential formula below.

3. The third formula also gets rid of the kinetic energy by computing the work.
Suppose the approximate force is related to a pseudo-potential $\tilde U$ as
$$
\tilde f = - \partial \tilde U/\partial x.
$$
The sum of the kinetic energy and the pseudo-potential is roughly conserved
by the integrator,
$$
\tilde U1 + K1 = \tilde U2 + K2.
$$
Then the change of the energy
$$
E2 - E1
 = (U2 + K2) - (U1 + K1)
 = (U2 - U1) - (\tilde U2 - \tilde U1)
 = (U2 - U1) + (1/2) (\tilde f2 + \tilde f1) \Delta x.
$$
This is referred to the work-formula below.


# Compilation

To compile the program, simply type
```
gcc dblwell.c -lm -o dblwell
```
or use `make`
```
make
```

# Running the program

To run the program, type
```
./dblwell
```
In the default case, HMC is not used, so there is no discontinuity at $xc$.
The distribution is saved in `hist.dat`.
This can be shown in by typing the following in Gnuplot
```
plot [-3:3] "hist.dat" u 1:(log($2)) w l
```

The discontinuity appears when we turn on HMC
```
./dblwell --hmc
```
Using the same the plotting command.

## Options

### Ways of computing the quantity $Q$

The default is the four-potential formula, or explicitly
```
./dblwell --hmc --hmcQ=U
```

To use the work formula:
```
./dblwell --hmc --hmcQ=W
```

To the total energy formula
```
./dblwell --hmc --hmcQ=E
```

### MD Integrators

The default MD integrator is velocity Verlet, or explicitly
```
./dblwell --hmc --int=vv
```

To the use the leapfrog formula, we set
```
./dblwell --hmc --int=leapfrog
```
Note that the leapfrog integrator,
when mixed with thermostat, is not time reversible.
So this algorithm is rather approximate.


### Multiple time step

An HMC step may be applied to multiple small MD steps.
```
./dblwell --hmc --block=20 --dt=0.4
```
In this case, the MD time step is $0.4/20 = 0.02$.

