#include "hmcrmsd/mdrun_main.h"

#include "gromacs/commandline/cmdlinemodulemanager.h"

int main(int argc, char *argv[])
{
    return gmx::CommandLineModuleManager::runAsMainCMain(argc, argv, &gmx_mdrun);
}
