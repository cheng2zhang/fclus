Overview
========


This directory contains a program `dblwell`
which is intended to show the application of
hybrid Monte Carlo (HMC), a.k.a. Hamiltonian Monte Carlo,
in sampling a distribution with discontinuity.

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

