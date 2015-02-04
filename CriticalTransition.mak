all:CriticalTransition.pdf

CriticalTransition.pdf: CriticalTransition.md margins.sty CriticalTransition.bib CriticalTransition.mak
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc -H margins.sty --bibliography CriticalTransition.bib --csl=oikos.csl CriticalTransition.md -o CriticalTransition.pdf 
	evince CriticalTransition.pdf		

CriticalTransition.docx: CriticalTransition.md margins.sty CriticalTransition.bib CriticalTransition.mak
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc -H margins.sty --bibliography CriticalTransition.bib --csl=oikos.csl CriticalTransition.md -o CriticalTransition.docx 
			
CriticalTransition_AMNAT.pdf: CriticalTransition_AMNAT.md margins.sty CriticalTransition.bib
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" CriticalTransition.bib
	pandoc -H margins.sty --bibliography CriticalTransition.bib --csl=the-american-naturalist.csl CriticalTransition_AMNAT.md -o CriticalTransition.pdf 
	pdftk CriticalTransition.pdf CriticalTransition_figures.pdf CriticalTransition_appendices.pdf output CriticalTransition_AMNAT.pdf
	evince CriticalTransition_AMNAT.pdf		
