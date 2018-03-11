work-rbras2018b: work-rbras2018b.Rnw rbras2018b.bib codes/*

	Rscript -e 'knitr::knit("work-rbras2018b.Rnw")'
	pdflatex work-rbras2018b.tex
	-bibtex work-rbras2018b.aux
	pdflatex work-rbras2018b.tex
	pdflatex work-rbras2018b.tex
	pdflatex abstract-pt-rbras2018b.tex
	pdflatex abstract-en-rbras2018b.tex
	make -s clean

simple:
	Rscript -e 'knitr::knit("work-rbras2018b.Rnw")'
	pdflatex work-rbras2018b.tex

abstract:
	pdflatex abstract-pt-rbras2018b.tex
	pdflatex abstract-en-rbras2018b.tex

clean:
	rm -f *.aux *.bbl *.blg *.brf *.idx *.ilg *.ind *.lof *.log \
	.*lol *.lot *.out *.toc *.synctex.gz
