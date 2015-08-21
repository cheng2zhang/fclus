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


Code prototype
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




