OPTS= -H margins.sty --bibliography ctnhm.bib --csl=oikos.csl 


%.pdf: %.md 
	pandoc  -H ctnhm_figures.sty -V geometry:margin=1.2cm  --latex-engine=xelatex $^ -o $@

all: ctnhm.pdf ctnhm_appendices.pdf SteadyStatePlots.pdf ctnhm.docx ctnhm_oikos_changes.pdf


SomeResults.pdf: SomeResults.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	pdftk ctnhm.pdf ctnhm_figures.pdf output ctnhm_figs.pdf
	evince SomeResults.pdf		

ctnhm.pdf: ctnhm.md 
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	evince ctnhm.pdf		

ctnhm.tex: ctnhm.md 
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	evince ctnhm.pdf		

ctnhm.docx: ctnhm.md 
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
			
ctnhm_AmNat.pdf: ctnhm_AmNat.md
	cp "/home/leonardo/BibTeX/Manuscritos-Critical Transition.bib" ctnhm.bib
	pandoc $< -o $@ $(OPTS)
	evince $@ 		

ctnhm_appendices.pdf: ctnhm_appendices.md
	pandoc -H ctnhm_appendices.sty ctnhm_appendices.md -o ctnhm_appendices.pdf 
	evince ctnhm_appendices.pdf		

ctnhm_oikos_changes.pdf: ctnhm.md

	git show 3bf59ed:ctnhm.md > oldms.md 
	pandoc oldms.md -o oldms.tex $(OPTS)
	pandoc ctnhm.md -o ctnhm.tex $(OPTS)
	latexdiff -t CTRADITIONAL oldms.tex ctnhm.tex > ctnhm_oikos_changes.tex
	pdflatex ctnhm_oikos_changes
	rm {ctnhm_oikos_changes,oldms,ctnhm}.tex
	rm oldms.md


