## Tables 

-------------------------------------------------------
  Side   No. Species   $\mu$   $\alpha$        $m$
                               (mean dist.)       
------- ------------- ------- -------------- ----------  
   100            64     0.2   2.04 (26.66)    0.00016    

   150                   	     2.08 (13.34)    0.0016     

   256                   	     2.18 (6.67)     0.016      

   512                                              
-------------------------------------------------------

Table: Parameters values used in the simulations of the neutral-hierarchical model. The parameter $\mu$ is the mortality rate; $\alpha$ is exponent of the inverse power law dispersal kernel, between brakets is the mean dispersal distance; and $m$ is the migration from the metacommunity. 

\newpage


-------------------------------------------------------
 Meta.       Mean Dist.       $m$      $p_c$   SE $p_c$    
 type 
----------- ----------- ---------  --------- ----------
    L             26.66  0.00016    0.001549 0.00003166

    L             26.66  0.0016     0.002865 0.00009039

    L             26.66  0.016      0.01244  0.00004243

    L             13.34  0.00016    0.001676 0.0002007 

    L             6.67   0.00016    0.001596 0.0002502 

    U             26.66  0.00016    0.001284 0.00007359

    U             26.66   0.0016    0.001778 0.0000229 

    U             26.66   0.016     0.01097  0.0006949 

    U             13.34  0.00016    0.001309 0.0001115 

    U             6.67   0.00016    0.001235 0.0002246 
-------------------------------------------------------

Table: Critical probabilities for infinite lattices. Where *Mean Dist.* is the mean dispersal distance,  $m$ is the migration parameter, $p_c$ the critical probability and SE $p_c$ the standar error of the critical probability.  


## Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SpanPvsRepl_T5000_64_side_meta.png}
\caption{Probability of Spanning cluster and critical probability for a spatial neutral/hierarchical model.  The columns represent two different metacommunity types: L, a metacommunity with logseries species abundance distribution (SAD); U, a metacommunity with a uniform SAD. The columns represent the side of simulation lattice.  The parameter $\rho$ is the control parameter representing the intensity of competition. The critical point was determined as the point where the spanning probability is 0.5 (vertical green line)}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/HvsRepl_T5000_64_side_meta_lin.png}
\caption{Shannon diversity index and critical probability for a spatial neutral/hierarchical model with different metacommunity types and different lattice size. Points are simulations of the model with the parameter $\rho$ representing the intensity of competition, the vertical green line is the critical point were the phase transition occurs}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/RichvsRepl_T5000_64_side_meta.png}
\caption{Richness and critical probability for simulated neutral/hierarchical communities as a function of the intensity of competition $\rho$. The columns represent different metacommunity types and the rows different grid *side*, where the total size of the grid is *side x side*. Points are repeated simulations of a the spatial model (n=30) and the vertical line is the critical point were the phase transition occurs}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SAD_T5000_64_512_meta.png}
\caption{Rank abundance diagrams for simulated neutral/hierarchical communities as a function of the intensity of competition $\rho$. Metacommunities have 64 species and two different species abundance distributions (SAD): logseries SAD; and uniform SAD. The size of the grid was 512*512 sites, the other parameters used are migration=0.00016, dispersal distance=26.66.}
\end{figure}

* Make figures of exponential rate from power law with exponential cutoff for patch distribution of species not spanning or most abuntant. Using $\rho$ near $p_c$ for that lattice size to show if there is a minimun.  

