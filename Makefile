default:
	$(MAKE) -C prog

clean:
	rm -rf $(bins) $(bins_d) a.out *.dat* gmon.out \
	  .*.un~ *~ */.*.un~ */*~ */*/.*.un~ */*/*~ \
	  r[0-9]*hs
	$(MAKE) -C prog clean
	rstrip.py -Rlv

excludes = --exclude=".*" --exclude="*~" --exclude="bak"

Dropbox: clean
	rsync -avzL $(excludes) * ~/Dropbox/fclus/

Bossman: clean
	rsync -avzL $(excludes) * /Bossman/cz1/fclus/

Bossman2: clean
	rsync -vzL $(excludes) * cz1@129.109.88.204:/Bossman/cz1/fclus/

