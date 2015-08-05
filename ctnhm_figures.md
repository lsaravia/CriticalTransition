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

-------------------------------------
 Metacomm.  Side  $\rho_c$  Critical
  Type                       Cluster 
---------- ------ -------- ----------
    L       100   0.00018   0.3093   

    L       150   0.00102   0.3777   

    L       256   0.00132   0.3521   

    L       512   0.0015    0.2862   

    U       100   0.00065   0.3564   

    U       150   0.00094   0.3609   

    U       256   0.0011    0.3583   

    U       512   0.00136   0.3448   
-------------------------------------

Table: Critical points $\rho_c$ and critical cluster size for the phase transition of a neutral/niche model. We used a range of simulation lattice sizes and two metacommunity types. The metacommunities  where L: logseries, U: uniform; *Side* was the side of the simulation lattice and the total size=$Side^2$. The other parameters were migration $m$=0.00016, mean dispersal distance=26.66. 

\newpage


---------------------------------------------------
 Metacomm.   Mean      $m$     $\rho_c$  SE$\rho_c$    
 type      Distance 
---------- -------- --------  --------- -----------
    L         26.66  0.00016   0.001549  0.00003166
 
    L         26.66  0.0016    0.002865  0.00009039
 
    L         26.66  0.016     0.01244   0.00004243
 
    L         13.34  0.00016   0.001676  0.0002007 
 
    L         6.67   0.00016   0.001596  0.0002502 
 
    U         26.66  0.00016   0.001284  0.00007359
 
    U         26.66   0.0016   0.001778  0.0000229 
 
    U         26.66   0.016    0.01097   0.0006949 
 
    U         13.34  0.00016   0.001309  0.0001115 
 
    U         6.67   0.00016   0.001235  0.0002246 
---------------------------------------------------

-----------------------------------------------------------------
Metacomm.  Metacomm.  Mean           $m$    $\rho_c$   SE$\rho_c$    
 species     type     Distance
--------- ---------- ---------- ---------  ---------- -----------
   64        L           6.67    0.00016    0.0016        0.00025
     
   64        L           13.34   0.00016    0.00168       0.0002 
     
   64        L           26.66   0.00016    0.00156       0.00004
     
   64        L           26.66    0.0016    0.00286       0.00009
     
   64        L           26.66   0.01596    0.01244       0.00004
     
   64        U           6.67    0.00016    0.00124       0.00022
     
   64        U           13.34   0.00016    0.00131       0.00011
     
   64        U           26.66   0.00016    0.00128       0.00007
     
   64        U           26.66    0.0016    0.00178       0.00002
     
   64        U           26.66   0.01596    0.01151       0.00102
     
   256       L           26.66   0.00016    0.00189       0.00015
     
   256       U           26.66   0.00016    0.00161       0.00013
-----------------------------------------------------------------

Table: Critical points $\rho_c$ for infinite lattices. Where *Mean Dist.* is the mean dispersal distance,  $m$ is the migration parameter, $p_c$ the critical probability and SE $p_c$ the standar error of the critical probability.  


## Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SpanPvsRepl_T5000_64_side_meta.png}
\caption{Probability of Spanning cluster for a spatial neutral/niche model as a function of the intensity of competition $\rho$. The columns represent two different metacommunity types: L, a metacommunity with logseries species abundance distribution (SAD); U, a metacommunity with a uniform SAD. The columns represent the side of the simulation lattice.  The vertical green line is the critical point,  the value for parameter $\rho$ where a phase transition between neutral and niche phases occurs. The critical point was determined as the point where the spanning probability is 0.5, the other parameters used were $m$=0.00016, dispersal distance = 26.66}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/HvsRepl_T5000_64_side_meta_lin.png}
\caption{Shannon diversity index and critical probability for a spatial neutral/niche model  as a function of the intensity of competition $\rho$. Columns represent with metacommunity types: Logseries is a metacommunity with logseries species abundance distribution (SAD), and the Uniform metacommunity have a uniform SAD, both with 64 species. Rows represent different lattice sizes. Points are independent simulations of the model. The parameter $\rho$ representing the intensity of competition is the control parameter and the vertical green line is the critical point were the phase transition occurs. Other parameters used were $m$=0.00016, dispersal distance = 26.66.}
\end{figure}



\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/RichvsRepl_T5000_64_side_meta.png}
\caption{Richness and critical probability for simulated neutral/niche communities as a function of the intensity of competition $\rho$. The columns represent different metacommunity types and the rows different grid *side*, where the total size of the grid is *side x side*. Points are repeated simulations of a the spatial model (n=30) and the vertical line is the critical point were the phase transition occurs}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SAD_T5000_64_512_meta.png}
\caption{Rank abundance diagrams for simulated neutral/hierarchical communities as a function of the intensity of competition $\rho$. Metacommunities have 64 species and two different species abundance distributions (SAD): logseries SAD; and uniform SAD. The size of the grid was 512*512 sites, the other parameters used are migration=0.00016, dispersal distance=26.66.}
\end{figure}

* Make figures of exponential rate from power law with exponential cutoff for patch distribution of species not spanning or most abuntant. Using $\rho$ near $p_c$ for that lattice size to show if there is a minimun.  

