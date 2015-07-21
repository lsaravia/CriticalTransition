OPTS= -H margins.sty --bibliography CriticalTransition.bib --csl=plos.csl 


%.pdf: %.md 
	pandoc  -V geometry:margin=1cm  --latex-engine=xelatex $^ -o $@

all: CriticalTransition.pdf SomeResults.pdf


SomeResults.pdf: SomeResults.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc $< -o $@ $(OPTS)
	evince SomeResults.pdf		

CriticalTransition.pdf: CriticalTransition.md makefile
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc $< -o $@ $(OPTS)
	evince CriticalTransition.pdf		
