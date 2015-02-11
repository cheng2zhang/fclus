prj=hmc

figs=

$(prj).pdf: $(prj).tex
	pdflatex $<
	pdflatex $<

clean:
	rm -rf *.log *.aux *~
	rstrip.py -Rlv

$(prj).zip: $(prj).tex $(prj)doc.pdf
	       Makefile \
	zip -r $@ $^

zip: $(prj).zip

