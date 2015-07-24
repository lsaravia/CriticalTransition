all:ctnhm.pdf

ctnhm.pdf: ctnhm.md margins.sty ctnhm.bib ctnhm.mak
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc -H margins.sty --bibliography ctnhm.bib --csl=oikos.csl ctnhm.md -o ctnhm.pdf 
	evince ctnhm.pdf		

ctnhm.docx: ctnhm.md margins.sty ctnhm.bib ctnhm.mak
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc -H margins.sty --bibliography ctnhm.bib --csl=oikos.csl ctnhm.md -o ctnhm.docx 
			
ctnhm_AMNAT.pdf: ctnhm_AMNAT.md margins.sty ctnhm.bib
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc -H margins.sty --bibliography ctnhm.bib --csl=the-american-naturalist.csl ctnhm_AMNAT.md -o ctnhm.pdf 
	pdftk ctnhm.pdf ctnhm_figures.pdf ctnhm_appendices.pdf output ctnhm_AMNAT.pdf
	evince ctnhm_AMNAT.pdf		
