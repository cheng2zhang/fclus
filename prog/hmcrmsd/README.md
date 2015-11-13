Overview
========

The directory contains a program `hmcrmsd`,
which is a molecular dynamics (MD) program
intended to sampling a flat histogram along
the root-mean-squared deviation (RMSD) from
a reference structure of a protein.

The program is a modification of the GROMACS
MD engine `mdrun`.  Thus, its use is similar.
The user should be able to compile GROMACS source.
Then, the files under this directory can be used to
supplement the GROMACS source code.


Installation
============

The source code under this directory is an independent module for GROMACS.

1. Install GROMACS 5.0
```
git clone git://git.gromacs.org/gromacs.git gromacs5.0
cd gromacs5.0
git checkout release-5-0
```

2. Copy this directory to `src/programs/`
3. Move `hmcrmsd_main.hpp` to `src/programs/hmcrmsd_main.cpp`. Note the change of the extension.
4. Append the following code to src/programs/CMakeLists.txt

```
file(GLOB MDRUN_SOURCES hmcrmsd/*.c hmcrmsd/*.cpp)
add_library(hmcrmsd_objlib OBJECT ${MDRUN_SOURCES})
add_executable(hmcrmsd $<TARGET_OBJECTS:hmcrmsd_objlib> hmcrmsd_main.cpp)
target_link_libraries(hmcrmsd libgromacs ${GMX_EXE_LINKER_FLAGS})
set(BINARY_NAME "hmcrmsd${GMX_BINARY_SUFFIX}")
set_target_properties(hmcrmsd PROPERTIES
    OUTPUT_NAME "${BINARY_NAME}"
    COMPILE_FLAGS "${OpenMP_C_FLAGS}")
install(TARGETS hmcrmsd DESTINATION ${BIN_INSTALL_DIR} COMPONENT hmcrmsd)
```

On stampede
--------------

```
cd work/gmx/gromacs5.0
mkdir buildicc && cd buildicc
module load intel/15.0.2 fftw3 cmake boost-mpi vmd
cmake .. -DCMAKE_CC_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
```

Note to load these modules when


System preparation
==================


1. Making the input PDB
-----------------------

Any protein PDB file can be used to as the input.

To make a PDB file for an ideal helix,
run the python `mkhelix.py` under this directory.
```
python mkhelix.py -n 12 --ter -o ala12.pdb
```

The option `-n 12` specifies the number of amino acid
residues. The default residue is alanine (ALA).

The option `--ter` adds a residue `ACE` on the N-terminal
at the beginning of the amino-acid chain
(an amino acid chain goes as `N-CA-CO-N-CA-CO-...`),
and a residue `NH2` on the C-terminal
at the end of the chain.
If only the N-terminal residue ACE is needed,
then change `--ter` to `--nter`;
if only the C-terminal residue NH2 is needed
then change `--ter` to `--cter`;
if the option is missing, no terminal residue is added.

The script is runnable, meanning that on a Linux machine,
the `python` at the beginning of the command can be dropped.



2. Using the script `simulpdb.py`
---------------------------------

This step converts the input PDB file to GROMACS files.

The python script `simulpdb.py` is a wrapper for GROMACS programs
to prepare initial files from a single PDB file.
The script is linked here under the subdirectory `simulpdb`

This step prepare the the following files for MD simulation in GROMACS.

File            | Description
----------------|------------------------
`init.gro`      | initial configuration
`topol.top`     | topology file
`mdrun.mdp`     | sample MD parameter file

Here is a set of typical commands of using the script
```
mkdir init
cd init
python path/to/simulpdb.py \
  --gmxexe=gromacs/build/root \
  -d 9
  --ff=amber03
  --sver=5.0
  my.pdb
cd ..
```

The argument `gromacs/build/root` for the option `--gmxexe=`
is `~/lwork/gmx/gromacs5.0/buildgcc` for the office computer.

The option `-d 9` means that the protein is separated from
any of periodic mirror image by at least 9 angstroms
(0.9 nm) in the simulation box.
The number can be adjusted, with a larger number means
more water molecules in the simulation box.

The force field in the above example is AMBER03, `--ff=amber03`.
If the option is missing the force field is AMBER99SB-ILDN.

The option `--ver=5.0` tells the script the GROMACS version.


In the office computer, the test system is prepared under
`~/lwork/gmx/gromacs5.0/buildgcc/hmctest_ala12`.
The command is
```
python ~/lwork/gmx/user/code/python/simulpdb/simulpdb.py \
  --gmxexe=~/lwork/gmx/gromacs5.0/buildgcc \
  -d 12 --ff=amber03 --sver=5.0 \
  ~/lwork/fclus/prog/hmcrmsd/ala12.pdb
```

Note `-d 9` creates a box of roughly 2000 water molecules,
`-d 12` creates a box of roughly 3000 water molecules.


### Notes

The original version of `simulpdb.py` is located under
```
WORKROOT/gmx/user/code/python/simulpdb/simulpdb.py
```
where `WORKROOT` is either `~/work` or `~/lwork`.


3. Preparation for the MD simulation
------------------------------------

### Copy files to the running directory

Link or copy the files `init.gro`, `topol.top` and `mdrun.mdp` from the previous step
to the running directory.

Also copy the reference PDB file, which is `ala12.pdb` in our example,
and the HMC parameter file, which is `ala12.cfg`, to the running directory.

### Change MD parameters in `md.mdp`

Rename `mdrun.mdp` to `md.mdp` and

 * Increase the number of steps by changing `nsteps=`.
 * Modify the number of steps for the interval of configuration output `nstxtc=`.
 * Modify the XTC group, `xtcgroup=Protein` if only protein coordinates are of interest.
 * Modify the number of steps of logging, `nstlog`
 * Modify the number of steps of energy output, `nstenergy`
 * Change the random number seed, `gen_seed`

```
gromacs/build/root/bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
```

For the alanine example in the office computer
```
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
```
or on stampede
```
~/work/gmx/gromacs5.0/buildicc/bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
```


### HMC parameters

The HMC parameters are controlled by the `.cfg` file.
See the special section below for details.



4. Run the MD simulation
------------------------

Running the program `hmcrmsd` is almost the same as running
the GROMACS MD engine `mdrun`, the only difference is that
we need to include the flag `-cfg myhmc.cfg` to specify the
file for HMC parameters.

```
gromacs/build/root/bin/hmcrmsd -cfg myhmc.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```

For the alanine example in the office computer
```
../bin/hmcrmsd -cfg ala12.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```
```
~/work/gmx/gromacs5.0/buildicc/bin/hmcrmsd -cfg ala12.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```


Parameters in `.cfg` file
=========================


The parameters of the program `hmcrmsd` are controlled by
the configuration file, such as `ala12.cfg`.
This is text file using a format similar to that of the GROMACS `.mdp` file.

Each line of the file specifies an option
```
key = value
```
The format for `key` is case-insentitive, and it allows some hyphenations.
Thus `RMSDmin` is equivalent to `RMSD-min` and `rmsdmin`.

A comment line starts with `#` or `;`.

Below is a list of common options.

## Most common options

The following are the most common options in the configuration file
for the program `hmcrmsd`

  * Change the name of the reference PDB,
    ```
    PDB = ala12.pdb
    ```

  * Check the target atom group for the RMSD bias,
    ```
    RMSD-group = heavy
    ```
    The option `heavy` (meaning all non-hydrogen atoms on the protein)
    can be changed to `CA` (alpha-carbon atoms) or `all` (all atoms on
    the protein).  The option `all` should be used with caution,
    because not all PDB files furnish coordinates for hydrogen atoms.
    Compared to `CA`, the option `heavy` is preferred because as more
    atoms are included in the RMSD group, the gentler the bias force
    would be.

  * Change the target RMSD range and bin size
    ```
    RMSD-min = 0.10
    RMSD-max = 0.70
    RMSD-del = 0.01
    ```
    Generally, `RMSD-max` should increase with the protein size.
    The idea is to cover the RMSD range for transition.

  * Change the target RMSD distribution to 1/RMSD^a,
    with exponent given by
    ```
    wr-exponent = 1.0
    ```
    The default exponent is 0.0, meaning a flat histogram.


## Advanced options

The following options are less important.  We usually don't need to
worry about them too much.  But in case something goes wrong, they
may have to look into them.

  * Wang-Landau parameters.
    The potential energy surface is continuously updated on-the-fly.
    The magnitude of updating is initially set to a value of `lnf0`,
    given by
    ```
    WL-lnf0 = 1e-3
    ```
    The parameter should be large enough to quickly achieve
    an initial flat distribution along RMSD.  This would give
    the simulator some hope for the simulation.
    However, it should not be too large as it can blow up the system.
    In the long run the parameter must be reduced to approach a more
    equilibrium-style simulation.
    The reduction of the updating magnitude of `lnf` is done stagewise.
    Particularly, in the initial stages, we will periodically
    check if the distribution along RMSD is flat enough, that is if
    the histogram flatness, measured by ratio of the standard deviation
    of the histogram heights to the average height,
    is less than the value specified by
    ```
    WL-flatness = 0.3
    ```
    If this is so, we reduce `lnf` by a factor of
    ```
    WL-frac = 0.5
    ```
    In the long run, the above strategy is no long effective,
    and we will simply use the formula of `lnf = C/t`
    to change the updating factor, where `t` is the number of
    steps per bin, and the value of `C` is given by
    ```
    Invt-C = 1
    ```

  * Valid range of the mean force,
    ```
    mfmin = -2000
    mfmax = 2000
    ```
    The minimal and maximal mean-force value along RMSD.
    The GROMACS unit (kJ/mol/nm) is used.

  * Valid range of the mean force at the left boundary,
    ```
    mflmin = -2000
    mflmax = -100
    ```
    These are the valid mean force range for `RMSD < RMSD-min`.
    In the above example, we only allow a negative mean force
    so that the RMSD of the system would be pushed to the desired
    range.

  * Similarly, we have parameters for the valid range the mean force
    at the right boundary
    ```
    mfhmin = 100
    mfhmax = 2000
    ```

## Input/output options

  * The computed potential of mean force (PMF) is saved in the file
    ```
    fnvrmsd = vrmsd.dat
    ```
    The first column is the RMSD, the second the PMF.
    The third column is the normalized histogram.
    The fourth column is the unnormalized histogram.
    To plot the PMF using GNUplot, use
    ```
    plot "vrmsd.dat" u 1:2 w l
    ```
    This file is saved every
    ```
    nstrep = 10000
    ```
    MD steps.

  * The time series of the RMSD is saved in the file
    ```
    fnlog = rmsd.log
    ```
    which is saved every
    ```
    nstlog = 100
    ```
    MD steps.

  * The frequency of printing command-line messages is controlled by
    ```
    nstchat = 1000
    ```
    which can be reduced when running on a supercomputer.


## Power user options

  * Feeling lucky to skip tests
    ```
    lucky = 1
    ```

## Testing/debugging options



The following options are meant to for debugging and future development.
This means that most of the options do NOT work.  So do not use them!

  * Enabling explicit HMC
    ```
    explicit-HMC = 1
    ```
    This option enables the explicit HMC algorithm.
    This means no RMSD bias force, and the flat distribution
    is entirely enforced by Metropolis acceptance criterion.
    However, if the RMSD drifts out of the desired range,
    the RMSD bias force is activated to pull it back into
    the desired range.

  * Enabling mean-force based histogram flattenning
    ```
    bias-mf = 1
    ```
    Instead of using the Wang-Landau based method, use a
    mean-force based method to a achieve a flat histogram.
    This method does not work, because
    1) No exact mean-force is available.
    2) Large mean-force fluctuation.
    3) Possible programming bugs.




Command line options
=========================

Command line options are mostly the same as those for GROMACS `mdrun`

  * `-cfg my.cfg`
    This option is only for `hmcrmsd`, it specifies the HMC-specific parameters

  * `-cpi md.cpt`
    This option, same as that for `mdrun`,
    is used to continue a simulation.
    When this option is used, the bias potential,
    updating magnitude, and RMSD histogram are inherited
    from the previous simulation.

  * `-ntmpi 2`
    This option specifies the number of MPI threads.
    These threads pretend to be MPI nodes, and they
    participate in domain decomposition.

  * `-ntomp 2`
    This option specifies the number of OMP threads,
    these threads do not participate in domain decomposition.

Source code
===========

Below is a brief note on the source code.

Code template
-------------

The source code is based on GROMACS 5.0.7.
Particularly, it was created by copying the directory
```
src/programs/mdrun
```
of the GROMACS source code tree.

However, `repl_ex.c`, `repl_ex.h`, `membed.h`
are deleted for simplicity.

The following options are also deleted for simplicity.

* rerun
* shellfc
* vsite
* pull
* FEP, SimTemp
* FAHCORE


Code analysis
-------------

### Review of integrators

#### leapfrog integrators

 v += (f/m) dt
 x += v dt

#### velocity-verlet integrators

 v += (f/m) (dt/2)      etrtVELOCITY2
 x += v dt              etrtPOSITIION
 v += (f/m) (dt/2)      etrtVELOCITY1

### `md.c`

This file contains the most important code.
The outline is shown below.
The symbol `[+]` denotes the modification
made for the program `hmcrmsd`.
The outline is made specifically for the leapfrog algorithm.

### for leapfrog


* MD loop starts on line 593

* `dd_partition_system()`, line 685

  o bMasterState is usually FALSE
  o usually go into branch 3, line 9531 in mdlib/domdec.c.


* `do_force()`, line 779

  o Defined in `sim_util.c`
  o `bStateChanged` is TRUE
    * `gmx_pme_send_coordinates()`, line 1738
    * `dd_move_x()`, line 1751
      o copy coordinates, `buf[n]`, line 705 in domdec.c
    * `ns()`, line 1828
    * `do_force_lowlevel()`, line 1925
      o gmx_pme_do();
    * `dd_move_f()`, line 1983
    * `pme_receive_force_ener()`, line 2057
      get force from PME nodes


* [+] `gmxgo_rmsd_force()`, line 788

  o collect the current coordinates
  o compute the force from the RMSD bias and add it to the total force `f`


* [+] `gmxgo_hmcpushxf()`, line 793

  o push the position and force.


* `do_md_trajectory_writing()`, line 998

  o The function call collects x, v, f so that
    the master node has the complete coordinates now.


* `update_tcouple()`, line 1190

  o update temperature coupling
    + every node does it (even nonmaster)
    + does it only `inputrec->nsttcouple` steps
  o defined in `gromacs/mdlib/update.c` as a wrapper
  o calls `vrescale_tcoupl()` defined in `gromacs/mdlib/coupling.c`
    +  `vrescale_tcoupl()` only sets `ekind->tcstat[i].lambda`
    +  it doesn't touch the actual velocity array.
  o for leapfrog, leave the velocity as is, will be done latter
        in `update_coords( etrtPOSITION )`


* `update_coords( etrtPOSITION )`, line 1222

  o defined in `gromacs/mdlib/update.c`
  o for leapfrog
    + update x, v for the local state
    + calls `do_update_md()` in the same file,
    + the normal branch starts from line 233
    + `v = v * ekind->tcstat[i].lambda + f * dt`
    + `xprime = x + v * dt`
    + `x` is usually not changed yet


* `update_constraints()`, line 1227

  o mainly `x = xprime` with modifications
  o defined in `gromacs/mdlib/update.c`
  o calls `constrain()`, defined in `gromacs/mdlib/constr.c`
    + coordinates are only partially communicated
    + change `x`; `xprime` is applied to `x`: `copy_rvec(upd->xp[i], state->x[i]);`
    + however, PBC is not applied! The molecule is not necessarily whole.


* [+] `gmxgo_hmcpushv()`, line 1299
  o push the velocity


* [+] `gmxgo_hmcselect()`, line 1438

  o decide whether to accept or reject the state `x`.


* `dd_collect_state()`, line 1525


* `dd_partition_system()`, line 1531

  o `bNeedRepartition`


* MD loop ends on line 1671



### For velocity Verlet

* [+] line 590
  o enforce T-coupling at every step.

* MD loop starts on line 593

* `dd_partition_system()`, line 685

  o bMasterState is usually FALSE
  o usually go into branch 3, line 9531 in mdlib/domdec.c.


* `do_force()`, line 779

  o Defined in `sim_util.c`
  o `bStateChanged` is TRUE
    * `gmx_pme_send_coordinates()`, line 1738
    * `dd_move_x()`, line 1751
      o copy coordinates, `buf[n]`, line 705 in domdec.c
    * `ns()`, line 1828
    * `do_force_lowlevel()`, line 1925
      o gmx_pme_do();
    * `dd_move_f()`, line 1983
    * `pme_receive_force_ener()`, line 2057
      get force from PME nodes


* [+] `gmxgo_rmsd_force()`, line 788

  o collect the current coordinates
  o compute the force from the RMSD bias and add it to the total force `f`


* [+] `gmxgo_hmcpushxf()`, line 793

  o push the position and force.


* `update_coords( etrtVELOCITY1 )`, line 831

  o second half VV step of the previous MD step
  o `update_coords( etrtVELOCITY1 )`
  o `do_update_vv_vel()`
  o v += (f/m) dt/2


* `do_md_trajectory_writing()`, line 998

  o The function call collects x, v, f so that
    the master node has the complete coordinates now.


* `update_tcouple()`, line 1113

  o update temperature coupling for the previous MD step.
  o for VV, update the velocity in place
  o [+] modified the step index by -1

* ################### The previous step completes now #################

* `update_tcouple()`, line 1190

  o update temperature coupling
    + every node does it (even nonmaster)
    + does it only `inputrec->nsttcouple` steps
  o defined in `gromacs/mdlib/update.c` as a wrapper
  o calls `vrescale_tcoupl()` defined in `gromacs/mdlib/coupling.c`
    +  `vrescale_tcoupl()` only sets `ekind->tcstat[i].lambda`
    +  it doesn't touch the actual velocity array.
  o for VV, update the velocity


* `update_coords( etrtVELOCITY2 )`, line 1199

  o first half VV step of this MD step
  o `update_coords( etrtVELOCITY2 )`
  o v += (f/m) dt/2


* `update_coords( etrtPOSITION )`, line 1222

  o defined in `gromacs/mdlib/update.c`
  o xprime = x + v dt
  o calls `do_update_vv_pos()`


* `update_constraints()`, line 1227

  o mainly `x = xprime` with modifications
  o defined in `gromacs/mdlib/update.c`
  o calls `constrain()`, defined in `gromacs/mdlib/constr.c`
    + coordinates are only partially communicated
    + change `x`; `xprime` is applied to `x`: `copy_rvec(upd->xp[i], state->x[i]);`
    + however, PBC is not applied! The molecule is not necessarily whole.


* [+] `gmxgo_hmcpushv()`, line 1299
  o push the velocity


* `update_coords(etrtPOSITION)`, line 1326

  o update position for VVAK


* [+] `gmxgo_hmcselect()`, line 1438

  o decide whether to accept or reject the state `x`.


* `dd_collect_state()`, line 1525


* `dd_partition_system()`, line 1531

  o `bNeedRepartition`


* MD loop ends on line 1671



### snippets

#### print `state->flags`

The following prints the flags
```
{
  int est;

  for ( est = 0; est < estNR; est++ )
    if ( state->flags & (1<<est) )
      fprintf(stderr, "%d: %s\n", est, est_names[est]);
}
```

### Personal notes for potential modifications

 *  set the default `nstglobalcomm` to 1 in `mdrun.cpp`?
 *  enforcing Andersen thermostat for explict-HMC runs?



