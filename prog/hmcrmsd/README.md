Overview
========

This module is intended to sampling a flat histogram along


Installation
============

1. Install GROMACS 5.0
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

5. To change running parameters, edit `hmcrmsd.cfg` and copy to the running directory.



System preparation
==================


1. Making the input PDB
-----------------------

Any protein PDB file can be used to as the input.

To make a PDB file for the ideal helix, run the mkhelix.py
```
python mkhelix.py --ter -n 12 -o ala12.pdb
```

If only the C-terminal residue NH2 is needed,
then change `--ter` to `--cter`.


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

```
mkdir init
cd init
path/to/simulpdb.py \
  --gmxexe=gromacs/build/root \
  -d 9
  --ff=amber03
  my.pdb
cd ..
```

The argument `gromacs/build/root` for the option `--gmxexe=`
is `~/lwork/gmx/gromacs5.0/buildgcc` for the office computer.

The option `-d 9` means that the protein is separated from
any of periodic mirror image by at least 9 angstroms
in the simulation box.
The number can be adjusted, with a larger number means
more water molecules in the simulation box.

The force field in the above example is AMBER03, `--ff=amber03`.
If the option is missing the force field is AMBER99SB-ILDN.


In the office computer, the test system is prepared under
`~/lwork/gmx/gromacs5.0/buildgcc/hmctest_ala12`.
The command is
```
~/lwork/gmx/user/code/python/simulpdb/simulpdb.py \
  --gmxexe=~/lwork/gmx/gromacs5.0/buildgcc \
  -d 9 --ff=amber03 \
  ~/lwork/fclus/prog/hmcrmsd/ala12.pdb
```


### Notes

The original version of `simulpdb.py` is located under
```
WORKROOT/gmx/user/code/python/simulpdb/simulpdb.py
```
where WORKROOT is either `~/work` or `~/lwork`.


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
 * Modify the number of steps of logging

```
gromacs/build/root/bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
```

For the alanine example in the office computer
```
../bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
```

### Change HMC parameters in `.cfg` file

The following are the most common options in the configuration file
for the program `hmcrmsd`

 * Change the name of the reference PDB,
   `PDB = ala12.pdb`.

 * Check the target atom group for the RMSD bias,
   `RMSD-group = heavy`.
   The option `heavy` (meaning all non-hydrogen atoms on the protein)
   can be changed to `CA` (alpha-carbon atoms) or `all` (all atoms on
   the protein).  The option `all` should be used with caution,
   because not all PDB files furnish coordinates for hydrogen atoms.
   Compared to `CA`, the option `heavy` is preferred because as more
   atoms are included in the RMSD group, the gentler the bias force
   would be.

 * Change the target RMSD range and bin size
    `RMSD-min = 0.10`
    `RMSD-max = 0.70`
    `RMSD-del = 0.01`


4. Run the MD simulation
------------------------

```
gromacs/build/root/bin/hmcrmsd -cfg myhmc.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```

For the alanine example in the office computer
```
../bin/hmcrmsd -cfg ala12.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```

Code template
===============

The code is based on GROMACS 5.0.7.
The source code is based on
```
src/programs/mdrun
```

`repl_ex.c`, `repl_ex.h`, `membed.h` are deleted for simplicity.

The following options are also deleted for simplicity.

* rerun
* shellfc
* vsite
* pull
* FEP, SimTemp
* FAHCORE


Modifications
=============

 *  set the default `nstglobalcomm` to 1 in `mdrun.cpp`



Code analysis
=============

md.c
----

* MD loop starts on line 582
* `dd_partition_system()`, line 674
  o bMasterState is usually FALSE
  o usually go into branch 3, line 9531 in mdlib/domdec.c.

* `do_force()`, line 768
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


* `do_md_trajectory_writing()`, line 977
  o The function call collects x, v, f so that
    the master node has the complete coordinates now.

* `update_tcouple()`, line 1168
  o update temperature coupling
    + every node does it (even nonmaster)
    + does it only `inputrec->nsttcouple` steps
  o defined in `gromacs/mdlib/update.c` as a wrapper
  o calls `vrescale_tcoupl()` defined in `gromacs/mdlib/coupling.c`
    +  `vrescale_tcoupl()` only sets `ekind->tcstat[i].lambda`
    +  it doesn't touch the actual velocity array.

* `update_coords()`, line 1200
  o update x, v for the local state
  o defined in `gromacs/mdlib/update.c`
  o calls `do_update_md()` in the same file,
    + the normal branch starts from line 233
    + `v = v * ekind->tcstat[i].lambda + f * dt`
    + `xprime = x + v * dt`
    + `x` is usually not changed yet

* `update_constraints()`, line 1205
  o defined in `gromacs/mdlib/update.c`
  o calls `constrain()`, defined in `gromacs/mdlib/constr.c`
    + coordinates are only partially communicated
    + `xprime` is applied to `x`: `copy_rvec(upd->xp[i], state->x[i]);`
    + however, PBC is not applied! The molecule is not necessarily whole.

* `dd_collect_state()`, line 1430

* `dd_partition_system()`, line 1436
  o `bNeedRepartition`

* MD loop ends on line 1576

snippets
--------

### print state->flags

The following prints the flags
```
{
  int est;

  for ( est = 0; est < estNR; est++ )
    if ( state->flags & (1<<est) )
      fprintf(stderr, "%d: %s\n", est, est_names[est]);
}
```
