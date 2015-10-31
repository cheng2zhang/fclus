all: stamp

stampaddr=oo1@stampede.tacc.utexas.edu

stamp:
	rsync -avz --exclude=*~ --exclude=.*.swp --exclude=.*.un~ \
	    *.c *.h *.cpp *.hpp *.py *.cfg *.pdb *.md \
	    $(stampaddr):work/gmx/gromacs5.0/src/programs/hmcrmsd/
