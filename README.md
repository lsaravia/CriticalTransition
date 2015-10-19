# Biodiversity collapse in a phase transition between neutral and niche communities

Markdown version of the manuscript and R source code for the data analysis. 

**ctnhm.md** the manuscript in markdown 

**ctnhm.bib** the bibliography

**ctnhm.pdf** the pdf version without figures

**ctnhm_figures.md** Manuscript figures

**ctnhm_appendices.md** Apendices


To reproduce the results you need to compile the following C++ program available at <https://github.com/lsaravia/Neutral>:

+ ipsNeutralCont: with power dispersal


## R Markdown files:

I usually run this chunk by chunk because there are long simulations, and save ".RData" to keep the results.
  
Neutral_simulCT.Rmd: Simulations to find critical values, and plots of results.


## R functions

R/Neutral_fun.r: collection of functions used for the analysis

R/powerlaw: code to fit powerlaws from <http://tuvalu.santafe.edu/~aaronc/powerlaws/> (by Cosma Shalizi)


# License

Unless explicitly stated in the source code or in a folder README the following applies:

    The MIT License (MIT)

    Copyright (c) 2015 Leonardo A. Saravia

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.



