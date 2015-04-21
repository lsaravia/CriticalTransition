OPTS= -H margins.sty --bibliography CriticalTransition.bib --csl=plos.csl 

all: SomeResults.pdf 


SomeResults.pdf: SomeResults.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc $< -o $@ $(OPTS)
	evince CriticalTransition.pdf		
