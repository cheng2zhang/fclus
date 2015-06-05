prj=hmc

figs=ljmix1
pdffigs=$(patsubst %,fig/%.pdf, $(figs))
epsfigs=$(patsubst %,fig/%.eps, $(figs))
figgps=$(patsubst %,fig/%.gp, $(figs))

latexpdfs = hmc.pdf short.pdf
latexsrcs = $(patsubst %.pdf,%.tex, $(latexpdfs))

bibfiles = liquid.bib simul.bib

all: $(latexpdfs)

$(latexpdfs): %.pdf: %.tex $(bibfiles) $(pdffigs) Makefile
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

$(pdffigs): $(epsfigs) $(figgps)
	make -C fig

clean:
	$(RM) -rf *.log *.aux *.blg *Notes.bib *.bbl *~
	rstrip.py -Rlv

$(prj).zip: $(latexsrcs) $(latexpdfs) \
  		$(epsfigs) $(pdffigs) $(figgps) \
		Makefile \
	zip -r $@ $^

zip: $(prj).zip

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/fclus/doc/
