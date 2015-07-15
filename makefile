OPTS= -H margins.sty --bibliography CriticalTransition.bib --csl=plos.csl 


%.pdf: %.md 
	pandoc  -V geometry:margin=1cm  --latex-engine=xelatex $^ -o $@

all: CSDCPoster.pdf

CSDCPoster.pdf: CSDCPoster.md makefile
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc $< -o $@  --bibliography CriticalTransition.bib



SomeResults.pdf: SomeResults.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc $< -o $@ $(OPTS)
	evince SomeResults.pdf		
