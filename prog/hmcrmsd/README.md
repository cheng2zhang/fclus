Overview
========

This module is intended to sampling a flat histogram along


Install
=======

1. Install GROMACS 5.0
2. Copy this directory to `src/programs/`
3. Move `hmcrmsd_main.cpp` to `src/programs/`
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

* `update_constraints()`, line 1205
  o defined in `gromacs/mdlib/update.c`
  o calls `constrain()`, defined in `gromacs/mdlib/constr.c`
    + coordinates are only partially communicated

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