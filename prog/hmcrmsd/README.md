Overview
========

This module is intended to sampling a flat histogram along


Install
=======

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

To make a PDB file for the ideal helix, run the mkhelix.py
```
python mkhelix.py --cter -n 20 -o helix20.pdb
```


2. Using the script `simulpdb.py`
---------------------------------

The python script `simulpdb.py` is a wrapper for GROMACS programs
to prepare initial files from a single PDB file.
The script is linked here under the subdirectory `simulpdb`

This step prepare the the following files for GROMACS simulation.

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



### Notes

The original version of `simulpdb.py` is located under
```
WORKROOT/gmx/user/code/python/simulpdb/simulpdb.py
```
where WORKROOT is either `~/work` or `~/lwork`.


3. Preparation for MD simulations
----------------------------------

```
gromacs/build/root/bin/gmx grompp -f md.mdp -c init.gro -o md.tpr
gromacs/build/root/bin/hmcrmsd -cfg 1LE1.cfg -deffnm md -v -ntmpi 1 -ntomp 2
```



Code template
===============

The code is based on GROMACS 5.0.5.
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
