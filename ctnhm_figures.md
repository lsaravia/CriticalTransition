## Tables 

-------------------------------------------------------
  Side   Metacomm.     $\mu$   $\alpha$        $m$
         No. Species           (mean dist.)       
------- ------------- ------- -------------- ----------  
   128           320     0.2   2.04 (26.6)    0.0001    

   192           64      	     2.08 (13.3)    0.001     

   256            32     	     2.02 (53.3)    0.01      

                                                
-------------------------------------------------------

Table: Parameters values used in the simulations of the neutral-hierarchical model. Side, is the size of the side of the simulation grid. The parameter $\mu$ is the mortality rate; $\alpha$ is exponent of the inverse power law dispersal kernel, between brackets is the mean dispersal distance; and $m$ is the migration from the metacommunity. The units of the simulation grid and dispersal are in meters to make them comparable with field values.

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


## Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/SpanPvsRepl_T5000_64_side_meta.png}
\caption{Probability of Spanning cluster for a spatial neutral/niche model as a function of the intensity of competition $\rho$. The columns represent two different metacommunity types: Logseries, a metacommunity with logseries species abundance distribution (SAD); Uniform, a metacommunity with a uniform SAD. The columns represent the side of the simulation lattice, the total size is $side^2$.  The vertical green line is the critical point,  the value for parameter $\rho$ where a phase transition between neutral and niche phases occurs. The critical point was determined as the point where the spanning probability is 0.5, the other parameters used were $m$=0.00016, dispersal distance = 26.66}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/HvsRepl_T5000_64_side_meta_lin.png}
\caption{Shannon diversity index and critical point for a spatial neutral/niche model  as a function of the intensity of competition $\rho$. Columns represent with metacommunity types: Logseries is a metacommunity with logseries species abundance distribution (SAD), and the Uniform metacommunity have a uniform SAD, both with 64 species. Rows represent different lattice sizes. Points are independent simulations of the model. The parameter $\rho$ representing the intensity of competition is the control parameter and the vertical green line is the critical point were the phase transition occurs. Other parameters used were $m$=0.00016, dispersal distance = 26.66.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/RichvsRepl_T5000_64_side_meta.png}
\caption{Richness and critical point for simulated neutral/niche model communities as a function of the intensity of competition $\rho$. The columns represent different metacommunity types and the rows different grid *side*, where the total size of the grid is *side x side*. Points are repeated simulations of a the spatial model (n=30) and the vertical line is the critical point were the phase transition occurs}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{figs/discRateOther_T5000_64_256_meta.png}
\caption{Exponential decay rate $\lambda$ for patch size distributions as a function of the intensity of competition $\rho$. We fitted a power law with exponential cutoff to patch size distribution of species that does not have the biggest patch (Other MaxPatch) or species that are not the spanning species (Other Spanning), the points are the value of $\lambda$ and the continuous lines are fitted median regressions. The dashed line connects the $\lambda$ median for each $\rho$.  
We made 10 simulations for each $\rho$, metacommunities have 64 species and two different species abundance distributions (SAD): *L*, logseries SAD; and *U*, uniform SAD; the critical point for logseries is 0.0015, for uniform metacommunities is 0.0014. The size of the grid was 256*256 sites and the other parameters used are migration=0.00016, dispersal distance=26.66.}
\end{figure}
