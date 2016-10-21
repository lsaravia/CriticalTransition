# Appendices 

## Appendix 1: Barro colorado Island forest plot discretization procedure

To analyze the patch distribution of BCI plot we first have to discretized the positions of the trees to fit them in a lattice. In each position only one individual of a particular species can be present, this is the same assumption that we made for the model that we use in the paper above. 

We have to choose a length scale to make the discretization, if we intend to fit all the individuals of all species in a different site the scale should be around 0.10 m for this plot, as the plot has 1000m x 500 m, that would result in an big lattice of 10000x5000 sites with a great proportion of empty places. This will result in a majority of isolated sites with almost no patches. If we use a bigger scale e.g. of 0.5 m, more than 1 individual of possibly different species get in some of the sites, in these cases we have to decide which one will occupy the site. We establish that the one with greater dbh, no matter the species, will be the one that occupies the site, thus we are favoring the more mature individuals.

In this process we have to find the scale that give us the maximal occupation of the lattice without loosing the species structure of the community. Then the criteria to stop enlarging the scale is that the species abundance distribution (SAD) of the discretized lattice has not to be statistically different from the original SAD. To test this we use the Anderson-Darling statistic with a randomization procedure using the R package kSamples [1], this statistic has been proved powerful to detect different kinds of communities [2]. Using this procedure we obtained a scale of 1 m, thus we used a lattice of 1000x500 sites. 

1. Scholz F, Zhu A (2015) kSamples: K-Sample Rank Tests and their Combinations. Available: http://cran.r-project.org/package=kSamples.

1. Saravia LA (2015) A new method to analyse species abundances in space using generalized dimensions. Methods Ecol Evol 6: 1298–1310. Available: http://doi.wiley.com/10.1111/2041-210X.12417.


## Appendix Tables 

----------------------------------------
               $\alpha$
  $m$      (mean dist.)  $\theta$    $I$  
--------- ------------- ---------  -----
    0.01    2.08 (13.3)     155.7  121.1

   0.001    2.08 (13.3)     113.8  97.3

   0.0001   2.08 (13.3)     83.2   78.2

    0.01    2.04 (26.6)     561.2  342.5

   0.001    2.04 (26.6)     410.3  275.2

   0.0001   2.04 (26.6)     300.0  221.1

    0.01    2.02 (53.3)    2026.0  969.7

   0.001    2.02 (53.3)    1482.0  779.2

   0.0001   2.02 (53.3)    1083.0  626.1
----------------------------------------

Table: Equivalence of spatially explicit parameters $m$ and $\alpha$ (Mean dispersal distance) with spatially implicit neutral parameters $\theta$ and $I$.  


---------------------------------------------------------------------------
 Side   MetaType    beforepc          $S_{max}$               $RS_{max}$ 
------ ---------- ------------- -------------------- ----------------------
 128   Logseries  $rho < rho_c$        0.938                 0.973         

 128   Logseries  $rho > rho_c$        0.147                 0.230         

 128    Uniform   $rho < rho_c$        0.951                 0.985         

 128    Uniform   $rho > rho_c$        0.127                 0.205         

 192   Logseries  $rho < rho_c$        0.934                 0.972         

 192   Logseries  $rho > rho_c$        0.020                 0.040         

 192    Uniform   $rho < rho_c$        0.942                 0.978         

 192    Uniform   $rho > rho_c$        0.031                 0.056         

 256   Logseries  $rho < rho_c$        0.928                 0.967         

 256   Logseries  $rho > rho_c$        0.002                 0.008         

 256    Uniform   $rho < rho_c$        0.940                 0.974         

 256    Uniform   $rho > rho_c$        0.002                 0.009         
---------------------------------------------------------------------------

Table: Size of the largest patch relative to the total area $S_{max}$ before the critical point $\rho < \rho_c$ and after the critical point $\rho > \rho_c$, and for the largest patch relative to the total species area $RS_{max}$. The parameters used were the specified in the first row of table 1. 





---------------------------------------------------------------------
Metacomm.  Metacomm.  Mean           $m$  $\rho_c^\infty$  SE$\rho_c$    
 species     type     Distance
--------- ---------- ---------- --------- --------------- -----------
   16         L          26.66    0.0001       0.00017      0.00003

   16         U          26.66    0.0001       0.00026      0.00001

   64         L          26.66    0.0001       0.00029      0.00001

   64         U          26.66    0.0001       0.00026      0.00000

   320        L          53.33    0.0001       0.00028      0.00002

   320        U          53.33    0.0001       0.00026      0.00001

   320        L          26.66    0.0001       0.00026      0.00002

   320        U          26.66    0.0001       0.00024      0.00001

   320        L          13.34    0.0001       0.00027      0.00000

   320        U          13.34    0.0001       0.00026      0.00001

   320        L          26.66    0.0010       0.00052      0.00008

   320        U          26.66    0.0010       0.00062      0.00007

   320        L          26.66    0.0100       0.00646      0.00000

   320        U          26.66    0.0100       0.00640      0.00000
---------------------------------------------------------------------

Table: Critical points $\rho_c^\infty$ for infinite lattices. Where *Mean Dist.* is the mean dispersal distance,  $m$ is the migration parameter, and SE $p_c$ the standard error of the critical point.  

\newpage

-----------------------------------------------------------
 Metacomm.   Variable     Delta          Delta    Relative 
 Type                   $\rho_c^\infty$ Variable Variation
---------- ------------ --------------- -------- ----------
    L       Dispersal      0.05            0.75    0.07    

    U       Dispersal      0.07            0.75    0.09    

    L        MetaNsp       0.42            0.95    0.44    

    U        MetaNsp       0.07            0.95    0.07    

    L          $m$         0.96            0.99    0.97    

    U          $m$         0.96            0.99    0.97    
-----------------------------------------------------------

Table: Relative variation of the critical point for infinite lattices $\rho_c^\infty$. The parameters are mean dispersal distance (Dispersal), metacommunity number of species (MetaNsp), and migration ($m$). See methods for details about simulations and table 1 for the ranges of parameters. We used metacommunities with two different species abundance distributions (SAD): *L* logseries SAD; and *U* uniform SAD.

\newpage

--------------------------------------
 model       type       n   Frequency 
------- -------------- --- -----------
NoModel     Spanning   121    0.25    

  Pow       Spanning   40     0.08    

  Pow       MaxPatch   65     0.13    

PowExp      Spanning   42     0.09    

PowExp      MaxPatch   221    0.46    

Exp         Spanning    0     0.00

Exp         MaxPatch    0     0.00
-------------------------------------

Table: Proportion of best models for patch size distributions from simulated neutral/niche model communities. We fitted 3 models to the patch distributions: exponential, power law (Pow) and power law with exponential cutoff (PowExp). The best model was selected using the Akaike information criteria. We made 30 simulations in a range of $\rho$ (see methods) and we used the following parameters: metacommunities have 64 species and two different species abundance distributions (SAD): logseries SAD; and uniform SAD; The size of the grid was 256*256 sites, migration=0.00016, and dispersal distance=26.66.

\newpage


## Appendix Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/{MaxClusRepl_T20000_320_side256_meta_m0.0001}.png}
\caption{Largest patch for a spatial neutral/niche model as a function of the intensity of competition $\rho$. The columns represent two different metacommunity types: Logseries, a metacommunity with logseries species abundance distribution (SAD); Uniform, a metacommunity with a uniform SAD. The rows represent the largest patch relative to total area $S_{max}$ and the largest patch relative to the species area $RS_{max}$. The vertical line is the critical point, which is the value for parameter $\rho$ where a phase transition between neutral and niche phases occurs. The parameters used were: side of the simulation lattice was 256 sites, the number of species in the metacomunity was 320, the metacommunity migration $m$=0.0001 and the dispersal distance = 26.66.}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SAD_T20000_320_256_meta_m0001.png}
\caption{Rank abundance diagrams (RADs) for simulated neutral/niche model communities as a function of the intensity of competition $\rho$. Except for $\rho=0$ the values in the legend are upper limits. The RADs are averages of 50 simulations. Metacommunities have 320 species and two different species abundance distributions (SAD): logseries SAD (L); and uniform SAD (U); the black line with $\rho=0.0002$ is the closest SAD previous to the critical point. The size of the grid was 256*256 sites and the other parameters used are  the migration rate $m$=0.0001, and dispersal distance=26.66.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/ClusBestModel_T20000_64_512_rho.png}
\caption{Proportion of best models as function of the competition intensity ($\rho$) for patch size distributions of the species with the largest patch or the spanning patch. We fitted fitted a power law (Pow), a power law with exponential cutoff, and an exponential models to distribution of patch sizes and selected the best model using the Akaike criterion. We made 30 simulations for each $\rho$ and metacommunity type. Metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 512*512 sites and the other parameters used were migration=0.0001, dispersal distance=26.66.  The critical point $\rho_c$ is 0.00029 for Logseries and 0.00026 for uniform metacommunities.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/AlfaVsRho_T20000_64_512_meta.png}
\caption{Power law exponent $\alpha$ for patch size distributions as a function of the intensity of competition $\rho$. We show the exponents of the two models selected by the Akaike criterion: the power law (Pow) and a power law with exponential cutoff. The patch size distribution corresponds to the species that has the largest patch or the species that percolate and form a spanning cluster. We made 30 simulations for each $\rho$ and metacommunity type. Metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 512*512 sites and the other parameters used were migration=0.0001, dispersal distance=26.66.  The critical point $\rho_c$ is 0.00029 for Logseries and 0.00026 for uniform metacommunities.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/ExpRateVsRho_T20000_64_512_meta.png}
\caption{Exponential decay rate of the Power law with exponential cutoff model for patch size distributions as a function of the intensity of competition $\rho$. We fitted a power law with exponential cutoff to patch size distribution of the species that has the biggest patch or the species that percolate and form a spanning cluster. We made 30 simulations for each $\rho$ and metacommunity type. Metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 512*512 sites and the other parameters used were migration=0.0001, dispersal distance=26.66.  The critical point $\rho_c$ is 0.00029 for Logseries and 0.00026 for uniform metacommunities.}
\end{figure}




