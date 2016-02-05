all: stamp

stampaddr=oo1@stampede.tacc.utexas.edu

stamp:
	rsync -avzL --exclude=*~ --exclude=.*.swp --exclude=.*.un~ \
	    ala12* *.mdp *.top *.gro \
	    $(stampaddr):scratch/hmcrmsd/ala12_99sb
