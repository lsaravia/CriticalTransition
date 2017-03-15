# Appendices 

Saravia, L.A., Momo, F.R. 2017 Biodiversity collapse and early warning indicators in a spatial phase transition between neutral and niche communities. Oikos.


## Appendix Tables 

-----------------------------------------------------
 Spatial     $\alpha$                     Non-spatial
  $m$      (mean dist.)   $\theta$    $I$     $m$ 
--------- -------------- ---------  ----- -----------
   0.01     2.08 (13.3)    155.7    121.1   0.0073   

  0.001     2.08 (13.3)    113.8    97.34   0.0059   

  0.0001    2.08 (13.3)    83.21    78.22   0.0048   

   0.01     2.04 (26.6)    561.2    342.5   0.0205   

  0.001     2.04 (26.6)    410.3    275.2   0.0165   

  0.0001    2.04 (26.6)     300     221.1   0.0133   

   0.01     2.02 (53.3)    2026     969.7   0.0559   

  0.001     2.02 (53.3)    1482     779.2   0.0454   

  0.0001    2.02 (53.3)    1083     626.1   0.0368   
----------------------------------------------------

Table: Equivalence of spatially explicit parameters $m$ and $\alpha$ (Mean dispersal distance) with spatially implicit neutral parameters $\theta$, $I$ and  non-spatial $m$.  We used the formulas 1a & 1b from Etienne & Rosindell (2011) that assume a fat tailed dispersal kernel and $128^2$ individuals. 

\newpage

--------------------------------------------------------------------------
 Side   MetaType   before/after       $S_{max}$               $RS_{max}$ 
------ ---------- --------------- ----------------- ----------------------
 128   Logseries  $\rho > \rho_c$        0.938                 0.973         

 128   Logseries  $\rho < \rho_c$        0.147                 0.230         

 128    Uniform   $\rho > \rho_c$        0.951                 0.985         

 128    Uniform   $\rho < \rho_c$        0.127                 0.205         

 192   Logseries  $\rho > \rho_c$        0.934                 0.972         

 192   Logseries  $\rho < \rho_c$        0.020                 0.040         

 192    Uniform   $\rho > \rho_c$        0.942                 0.978         

 192    Uniform   $\rho < \rho_c$        0.031                 0.056         

 256   Logseries  $\rho > \rho_c$        0.928                 0.967         

 256   Logseries  $\rho < \rho_c$        0.002                 0.008         

 256    Uniform   $\rho > \rho_c$        0.940                 0.974         

 256    Uniform   $\rho < \rho_c$        0.002                 0.009         
---------------------------------------------------------------------------

Table: Size of the largest patch relative to the total area $S_{max}$ before the critical point $\rho < \rho_c$ and after the critical point $\rho > \rho_c$, and for the largest patch relative to the total species area $RS_{max}$. The parameters used were the specified in the first row of table 1. 


\newpage


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

Table: Critical points $\rho_c^\infty$ for infinite lattices. Where the first column is the number of species in the metacommunity, the second  is the species abundance distribution of the metacommunity: *L* logseries and *U* uniform. *Mean Distance* is the mean dispersal distance,  $m$ is the migration parameter, and SE $p_c$ the standard error of the critical point.  


\newpage

--------------------------------------
 model       type       n   Frequency 
------- -------------- --- -----------
NoModel     Spanning   121    0.25    

  Pow       Spanning   40     0.08    

  Pow       MaxPatch   65     0.13    

PowExp      Spanning   42     0.08    

PowExp      MaxPatch   221    0.46    

Exp         Spanning    0     0.00

Exp         MaxPatch    0     0.00
-------------------------------------

Table: Proportion of best models for patch size distributions from simulated neutral/niche model communities. We fitted the distribution of patches of the species with the largest patch (MaxPatch) or the spanning patch (Spanning). We used 3 models: exponential (Exp), power law (Pow) and power law with exponential cutoff (PowExp), when there is not enough number of patches we did not fit any model (NoModel). The best model was selected using the Akaike information criteria. We made 30 simulations in a range of $\rho$ (see methods) and we used the following parameters: metacommunities have 64 species and two different species abundance distributions: logseries and uniform SAD. The size of the grid was 512*512 sites, migration=0.0001, and dispersal distance=26.66. 

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
\caption{Power law exponent $\alpha$ for patch size distributions as a function of the intensity of competition $\rho$. We show the exponents of the two models selected by the Akaike criterion: the power law (Pow) and a power law with exponential cutoff. The continuous lines unite the medians for each $\rho$. The patch size distribution corresponds to the species that has the largest patch or the species that percolate and form a spanning cluster. We made 30 simulations for each $\rho$ and metacommunity type. Metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 512*512 sites and the other parameters used were migration=0.0001, dispersal distance=26.66.  The critical point $\rho_c$ is 0.00029 for Logseries and 0.00026 for uniform metacommunities.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/ExpRateVsRho_T20000_64_512_meta.png}
\caption{Exponential decay rate of the Power law with exponential cutoff model for patch size distributions as a function of the intensity of competition $\rho$. The continuous line unites the medians for each $\rho$. We fitted a power law with exponential cutoff to patch size distribution of the species that has the biggest patch or the species that percolate and form a spanning cluster. We made 30 simulations for each $\rho$ and metacommunity type. Metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD. The size of the grid was 512*512 sites and the other parameters used were migration=0.0001, dispersal distance=26.66.  The critical point $\rho_c$ is 0.00029 for Logseries and 0.00026 for uniform metacommunities.}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/VarDeltaTRSmaxvsCP_T20000_320_256_meta.png}
\caption{Variance of temporal fluctuations of the largest patch species relative to the total abundance of the same species  $\Delta RS_{max}$. We simulated communities in the same time span than the simulations to determine the critical point---typically around 20000 time steps---we take the last 5000 and measure the patch sizes each 100 time steps. The communities that did not have a spanning patch were classified as "Before" the critical point, with a range of $\rho: 0 - 0.0004$. The communities that present a spanning patch in all the times are measured as "After" the critical point, with $\rho: 0.0004 - 1$. The communities where the spanning patch appears and disappears were classified as "Near" the critical point, with  $\rho: 0.0002 - 0.0004$. 
We made 10 simulations for each $\rho$ and two metacommunity types: "Logseries" species abundance distribution (SAD) and "Uniform" SAD. Metacommunities have 320 species, the size of the grid was 256*256 sites, migration from metacommunity was 0.0001, dispersal distance=26.66.}  
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SkewDeltaTRSmaxvsCP_T20000_320_256_meta.png}
\caption{Skewness of the temporal fluctuations of the largest patch species relative to the total abundance of the same species. We simulated communities in the same time span than the simulations to determine the critical point---typically around 20000 time steps---we take the last 5000 and measure the patch sizes each 100 time steps. The communities that did not have a spanning patch were classified as "Before" the critical point, whit a range of $\rho: 0 - 0.0004$. The communities that present a spanning patch in all the times are measured as "After" the critical point, with $\rho: 0.0004 - 1$. The communities where the spanning patch appears and disappears were classified as "Near" the critical point, with  $\rho: 0.0002 - 0.0004$. 
We made 10 simulations for each $\rho$ and two metacommunity types: "Logseries" species abundance distribution (SAD) and "Uniform" SAD. Metacommunities have 320 species, the size of the grid was 256*256 sites, migration from metacommunity was 0.0001, dispersal distance=26.66.}  
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discAlfaYear_BCI.png}
\caption{Power law exponent $\alpha$ for patch size distributions of the Barro Colorado Island forest plot as a function of the census year. The model fitted was a power law with exponential cutoff. The continuous line is a median regression (Slope: -0.0375, SE: 0.0145, t-value: -2.583, p-value: 0.049)}  
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discLamdaYear_BCI.png}
\caption{Exponential decay rate $\lambda$ for patch size distributions of the Barro Colorado Island forest plot as a function of the census year. The model fitted was a power law with exponential cutoff. The continuous line is a median regression (Slope: -0.0074, SE: 0.0033, t-value: -2.216, p-value: 0.078)}  
\end{figure}