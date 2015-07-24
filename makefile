OPTS= -H margins.sty --bibliography ctnhm.bib --csl=plos.csl 


%.pdf: %.md 
	pandoc  -V geometry:margin=1cm  --latex-engine=xelatex $^ -o $@

all: ctnhm.pdf ctnhm_figures.pdf


SomeResults.pdf: SomeResults.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	evince SomeResults.pdf		

ctnhm.pdf: ctnhm.md makefile
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	evince ctnhm.pdf		
