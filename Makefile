default:
	$(MAKE) -C prog

clean:
	$(MAKE) -C prog $@
	rm -rf .*.un~ *~ */.*.un~ */*~ */*/.*.un~ */*/*~ \
	  */*/*/.*.un~ */*/*/*~ \
	  r[0-9]*hs
	rstrip.py -Rlv

excludes = --exclude=".*" --exclude="*~" --exclude="bak" \
	   --exclude="*.dat" --exclude="*.pos" \
	   --exclude="tmp*" --exclude="_*"

Dropbox: clean
	rsync -avzL $(excludes) * ~/Dropbox/fclus/

Bossman: clean
	rsync -avzL $(excludes) * /Bossman/cz1/fclus/

Bossman2: clean
	rsync -vzL $(excludes) * cz1@129.109.88.204:/Bossman/cz1/fclus/

