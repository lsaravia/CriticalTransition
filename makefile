# makefile

%.html: %.md 
	pandoc --self-contained -s -S -i -t slidy -V slidy-url=slidy $^ -o $@
#	pandoc --self-contained $^ -o $@

%.pdf: %.md 
#	pandoc  -V geometry:margin=.5in  --latex-engine=xelatex $^ -o $@
	pandoc  -t beamer $^ -H Talks/Solabima_Charla.sty -V theme:Singapore --latex-engine=xelatex -o $@
	evince $@


all: Tandar_Talk.pdf
	
