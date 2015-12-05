# Appendices 

## Appendix 1: Barro colorado Island forest plot discretization procedure

To analyze the patch distribution of BCI plot we first have to discretized the positions of the trees to fit them in a lattice. In each position only one individual of a particular species can be present, this is the same assumption that we made for the model that we use in the paper above. 

We have to choose a length scale to make the discretization, if we intend to fit all the individuals of all species in a different site the scale should be around 0.10 m for this plot, as the plot has 1000m x 500 m, that would result in an big lattice of 10000x5000 sites with a great proportion of empty places. This will result in a majority of isolated sites with almost no patches. If we use a bigger scale e.g. of 0.5 m, more than 1 individual of possibly different species get in some of the sites, in these cases we have to decide which one will occupy the site. We establish that the one with greater dbh, no matter the species, will be the one that occupies the site, thus we are favoring the more mature individuals.

In this process we have to find the scale that give us the maximal occupation of the lattice without loosing the species structure of the community. Then the criteria to stop enlarging the scale is that the species abundance distribution (SAD) of the discretized lattice has not to be statistically different from the original SAD. To test this we use the Anderson-Darling statistic with a randomization procedure using the R package kSamples [1], this statistic has been proved powerful to detect different kinds of communities [2]. Using this procedure we obtained a scale of 1 m, thus we used a lattice of 1000x500 sites. 

1. Scholz F, Zhu A (2015) kSamples: K-Sample Rank Tests and their Combinations. Available: http://cran.r-project.org/package=kSamples.

1. Saravia LA (2015) A new method to analyse species abundances in space using generalized dimensions. Methods Ecol Evol 6: 1298–1310. Available: http://doi.wiley.com/10.1111/2041-210X.12417.


## Appendix Tables 


---------------------------------------------------------------------
Metacomm.  Metacomm.  Mean           $m$  $\rho_c^\infty$  SE$\rho_c$    
 species     type     Distance
--------- ---------- ---------- --------- --------------- -----------
   16         L          26.66    0.00016        0.00138     0.00012
   
   16         U          26.66    0.00016        0.00101     0.00004
   
   64         L          26.66    0.00016        0.00156     0.00004
   
   64         U          26.66    0.00016        0.00128     0.00007
   
   256        L          26.66    0.00016        0.00189     0.00015
   
   256        U          26.66    0.00016        0.00161     0.00013
   
   64         L          13.34    0.00016        0.00168     0.0002 
   
   64         U          13.34    0.00016        0.00131     0.00011
   
   64         L          6.67     0.00016        0.0016      0.00025
   
   64         U          6.67     0.00016        0.00124     0.00022
   
   64         L          26.66     0.0016        0.00286     0.00009
   
   64         U          26.66     0.0016        0.00178     0.00002
   
   64         L          26.66    0.01596        0.01244     0.00004
   
   64         U          26.66    0.01596        0.01151     0.00102
--------------------------------------------------------------------

Table: Critical points $\rho_c^\infty$ for infinite lattices. Where *Mean Dist.* is the mean dispersal distance,  $m$ is the migration parameter, and SE $p_c$ the standar error of the critical point.  

\newpage

-----------------------------------------------------------
 Metacomm.   Variable     Delta          Delta    Relative 
 Type                   $\rho_c^\infty$ Variable Variation
---------- ------------ --------------- -------- ----------
    L       Dispersal      0.07            0.75      0.1   

    U       Dispersal      0.06            0.75      0.08  

    L        MetaNsp       0.27            0.94      0.29  

    U        MetaNsp       0.37            0.94      0.39  

    L      Colonization    0.87            0.99      0.88  

    U      Colonization    0.89            0.99      0.9   
-----------------------------------------------------------

Table: Relative variation of the critical point for infinite lattices $\rho_c^\infty$. See methods for details about simulations and table 1 for the ranges of parameters. We used metacommunitis with two different species abundance distributions (SAD): *L* logseries SAD; and *U* uniform SAD.

\newpage

--------------------------------------
 model       type       n   Frequency 
------- -------------- --- -----------
NoModel    Spanning    91     0.22    

  Pow      Spanning    21     0.05    

  Pow      MaxPatch    13     0.03    

PowExp     Spanning    42      0.1    

PowExp     MaxPatch    137    0.33    

PowExp  Other MaxPatch 47     0.11    

PowExp  Other Spanning 60     0.15    
--------------------------------------

Table: Proportion of best models for patch size distributions from simulated neutral/niche model communities. We fitted 3 models to the patch distributions: exponential, power law (Pow) and power law with exponential cutoff (PowExp). The best model was selected using the Akaike information criteria. We made 10 simulations in a range of $\rho$ (see methods) and we used the following parameters: metacommunities have 64 species and two different species abundance distributions (SAD): logseries SAD; and uniform SAD; The size of the grid was 256*256 sites, migration=0.00016, and dispersal distance=26.66.

\newpage

-------------------------------------------------------------------------------
 Metacomm.  param    species       Tau   Value   Std. Error   t value   pvalue 
 Type
---------- ------ --------------  ----- ------- ------------ --------- --------
    L      lambda Other MaxPatch  0.25  -35.56     15.88       -2.24     0.03  

    L                              0.5  -31.36     12.97       -2.42     0.02  

    L                             0.75  -14.07      8.9        -1.58     0.12  

    U                             0.25   -72.9     24.63       -2.96     0.01  

    U                              0.5  -72.14      8.56       -8.42      0    

    U                             0.75  -71.44     15.39       -4.64      0    

    L             Other Spanning  0.25   25.34      15.1       1.68      0.1   

    L                              0.5   27.48     12.44       2.21      0.03  

    L                             0.75   35.46     12.45       2.85      0.01  

    U                             0.25   30.89      4.51       6.86       0    

    U                              0.5   23.52      6.9        3.41       0    

    U                             0.75   21.77      9.04       2.41      0.02  

    L       alpha Other MaxPatch  0.25   0.96      18.64       0.05      0.96  

    L                              0.5   -5.6      34.99       -0.16     0.87  

    L                             0.75   53.34      26.1       2.04      0.05  

    U                             0.25   15.94     21.22       0.75      0.46  

    U                              0.5   61.92     28.37       2.18      0.04  

    U                             0.75   85.21     63.78       1.34      0.19  

    L             Other Spanning  0.25  -28.36     16.91       -1.68     0.1   

    L                              0.5  -28.58      13.2       -2.17     0.04  

    L                             0.75  -20.06     13.44       -1.49     0.14  

    U                             0.25  -25.49     12.39       -2.06     0.04  

    U                              0.5  -23.13      10.6       -2.18     0.03  

    U                             0.75  -19.96      7.92       -2.52     0.01  

    L       alpha    MaxPatch     0.25   17.55      4.01       4.37       0    

    L                              0.5   15.62      8.75       1.79      0.08  

    L                             0.75   34.03     11.15       3.05       0    

    U                             0.25   3.26      10.39       0.31      0.75  

    U                              0.5   15.79     13.29       1.19      0.24  

    U                             0.75   36.25     26.98       1.34      0.18  

    L                Spanning     0.25   9.07      40.38       0.22      0.82  

    L                              0.5   36.98     67.01       0.55      0.58  

    L                             0.75   151.7      97.7       1.55      0.13  

    U                             0.25   151.4     47.05       3.22       0    

    U                              0.5   214.3      86.5       2.48      0.02  

    U                             0.75   145.5     213.2       0.68      0.5   

    L      lambda    MaxPatch     0.25   -8.73      3.68       -2.37     0.02  

    L                              0.5  -12.58      5.21       -2.41     0.02  

    L                             0.75   -6.11      7.65       -0.8      0.43  

    U                             0.25  -67.67      7.76       -8.72      0    

    U                              0.5  -75.12     17.14       -4.38      0    

    U                             0.75  -64.21      7.85       -8.18      0    

    L                Spanning     0.25   0.17       1.77        0.1      0.92  

    L                              0.5   -1.41      5.17       -0.27     0.79  

    L                             0.75   -1.41     11.11       -0.13     0.9   

    U                             0.25   2.92       2.78       1.05      0.3   

    U                              0.5   8.36       6.99        1.2      0.24  

    U                             0.75   6.66      29.32       0.23      0.82  
-------------------------------------------------------------------------------

Table: Quantile regression of patch distribution model parameters vs. $\rho$, the intensity of competition.  We fitted a power law with exponential cutoff to patch size distribution which has two parameters: *alpha* is the power exponent, and *lambda* is the exponential decay rate (See methods for functional formulas). We fitted the model to patches of the species that has the biggest patch (MaxPatch), species that form a spanning cluster (Spanning), all the species that are not the MaxPatch (Other MaxPatch) and all species that are not the spanning species (Other Spanning). We fitted 3 quantiles *Tau*=0.25,0.50 and 0.75 and we used a bootstraping procedure to assess significance. We made 10 simulations for each $\rho$, metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 256*256 sites and the other parameters used were migration=0.00016, dispersal distance=26.66.


--------------------------------------------------------------------
     type       param   Tau   Value   Std. Error   t value   pvalue 
-------------- ------- ----- ------- ------------ --------- --------
   MaxPatch     alpha   0.5     0        0.04       0.03      0.98  

   MaxPatch    lambda   0.5     0        0.01       -0.21     0.84  

Other MaxPatch lambda   0.5     0         0           0        1    

Other MaxPatch  alpha   0.5   0.01       0.01       1.25      0.27  
--------------------------------------------------------------------

Table: Median regression of patch distribution model parameters vs. Year for BCI data.  We fitted a power law with exponential cutoff to patch size distribution which has two parameters: *alpha* is the power exponent, and *lambda* is the exponential decay rate (See methods for functional formulas). We fitted the model to patches of the species that has the biggest patch (MaxPatch), and all the species that are not the MaxPatch (Other MaxPatch). We used a bootstraping procedure to assess significance. 


## Appendix Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SAD_T5000_64_512_meta.png}
\caption{Rank abundance diagrams (RADs) for simulated neutral/niche model communities as a function of the intensity of competition $\rho$. Except for $\rho=0$ the values in the legend are upper limits. The RADs are averages of 30 simulations. Metacommunities have 64 species and two different species abundance distributions (SAD): logseries SAD (L); and uniform SAD (U); the critical point for logseries is 0.0013 and 0.0011 for uniform metacommunities. The size of the grid was 256*256 sites and the other parameters used are migration=0.00016, dispersal distance=26.66.}
\end{figure}


\begin{figure}[H]
\centering
\begin{tabular}{cc}
\subfloat{\includegraphics[width=3.3in]{figs/CritProb_Migr_T5000_64_meta.png}} &
\subfloat{\includegraphics[width=3.3in]{figs/CritProb_Disp_T5000_64_meta.png}} \\
\subfloat{\includegraphics[width=3.3in]{figs/CritProb_MetaNsp_T5000_64_meta.png}} &
\end{tabular}
\caption{Critical point for infinite lattices $\rho_c^\infty$ in function of migration from metacommunity,local dispersal distance and number of species in the metacommunity. The parameters used in the simulations are specified in table 1. Two kinds of metacommunities with different species abundance distributions (SAD) were used: *L*, logseries SAD; and *U*, uniform SAD.}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discExpSpan_T5000_64_256_meta.png}
\caption{Power law exponent $\alpha$ for patch size distributions as a function of the intensity of competition $\rho$. We fitted a power law with exponential cutoff to patch size distribution of species that has the biggest patch (MaxPatch) or species that percolate and form a  spanning cluster (Spanning). We made 10 simulations for each $\rho$, metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD; the critical point for logseries is 0.0013, for uniform metacommunities is 0.0011. The size of the grid was 256*256 sites and the other parameters used were migration=0.00016, dispersal distance=26.66.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discRateSpan_T5000_64_256_meta.png}
\caption{Exponential decay rate $\lambda$ for patch size distributions as a function of the intensity of competition $\rho$. We fitted a power law with exponential cutoff to patch size distribution of species that has the biggest patch (MaxPatch) or species that percolate and form a  spanning cluster (Spanning). We made 10 simulations for each $\rho$, metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD; the critical point for logseries is 0.0013, for uniform metacommunities is 0.0011. The size of the grid was 256*256 sites and the other parameters used were migration=0.00016, dispersal distance=26.66.}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discExpOther_T5000_64_256_meta.png}
\caption{Power law exponent $\alpha$ for patch size distributions as a function of the intensity of competition $\rho$. We fitted a power law with exponential cutoff to patch size distribution of species that does not have the biggest patch (Other MaxPatch) or species that are not the spanning species (Other Spanning). We made 10 simulations for each $\rho$, metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD; the critical point for logseries is 0.0013, for uniform metacommunities is 0.0011. The size of the grid was 256*256 sites and the other parameters used were migration=0.00016, dispersal distance=26.66.}
\end{figure}



