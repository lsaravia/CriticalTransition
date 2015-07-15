
# The critical transition between hierarchical competition and neutrality in ecological communities

# A spatial phase transition between neutral and niche ecological communities

**Leonardo A. Saravia**, Ph.D.

Instituto de Ciencias Básicas

Universidad Nacional de General Sarmiento

J.M. Gutierrez 1159 (1613), Los Polvorines

Buenos Aires, Argentina.

<lsaravia@ungs.edu.ar>



## Abstract

Ecological communities are shaped by a mixture of niche-based and neutral ecological processes. In neutral communities individuals are identical except for is abundance in a metacommunity and the dynamics is dominated by demographic stochasticity and dispersal. Niche communities can be represented by competitive hierarchies were replacement of competitive superior species occurs. Here I present the most simple model representing the continuum of neutral-niche communities. I add to a classical neutral model a competitive hierarchy where each species have a probability $\rho$ of replace another species that measures competitive intensity. 
In this model shows a phase transition between neutral an niche communities that produce a spanning cluster of the most competitive species is observed. The diversity (measured with Shannon index) have and abrupt fall but richness shows a slower decrease. The critical point is observed at a very low value of the competitive intensity, thus weak competition can have a strong effect on communities. Some recent studies suggest that most communities present weak competitive interactions thus they may be near the critical point, understanding the dynamics this phase transition will be very important to prevent shifts in community structure.  

## Abstract (to present the model and this lack of trade-off)

Ecological communities are shaped by a mixture of niche-based and neutral ecological processes. In neutral communities individuals are identical except for is abundance in a metacommunity and the dynamics is dominated by demographic stochasticity and dispersal. Niche communities can be represented by competitive hierarchies and to maintain diversity some form of trade-off is assumed. Here I present the most simple model representing the continuum of neutral-niche communities. I add to a classical neutral model a competitive hierarchy where each species have a probability $\rho$ of replace another species. 
The model shows a phase transition between neutral an niche communities that produce a biodiversity collapse. This phase transition is characteristic of some ecological systems and has been termed robust criticality. Here we describe for the first time the observation of this kind of criticality for a community assembly model.  

## Abstract

The dynamics of ecological communities can be described by two contrasting models: the first assumes that the individuals of all species are identical and do not have competitive interactions, the second assumes that species are different, adapted to particular habitat conditions and have strong interactions. These represent extremes of a continuum: the first is the neutral and the later the niche model of communities. Real communities are actually a mixture of both dynamics. Here I study the simplest model of neutral-niche communities where niche dynamics is represented as a competitive hierarchy. The competition intensity is represented as a parameter that modulates the transition between these extremes. I use a stochastic cellular automata to show that there is a phase transition between the neutral and niche model with a spanning cluster formed by the most abundant species. The transition implies a sharp fall of species diversity (Shannon entropy) but the richness (number of species) shows a gentle decline with increasing competitive intensity. The critical point is at a very low value of competitive intensity and the same was reported for different real communities. Thus it will be possible that a mechanism of self-organization is driving communities to criticality, that mechanism could be similar but not equal to self-organizing criticality.


## Introduction

Much effort has been devoted to understand the mechanisms of community assembly and dynamics. In principle the emphasis were on deterministic mechanisms based on niche differences between species; the niche theory assumes that different species are regulated by different environmental factors and infer that diversity originates from spatial and temporal environmental heterogeneity [@Hutchinson1957; @Tilman1982; @Chesson2000]. More recently the emphasis shifted to stochastic mechanisms in the form of the Neutral theory of Biodiversity and Biogeography [@Hubbell2001]. The neutral theory assumes that individuals of all species are equivalent and it proposes that diversity originates from a balance between immigration, speciation, and extinction. Neutral theory has been proposed as a parsimonious formulation that can provide new insight into the patterns of community assembly [@Hubbell2005], besides this simplification it can predict some community metrics very well [@Volkov2007; @Rosindell2012], mainly the species abundance distribution (SAD).

Finally arises a unified view that accepts that both kinds of mechanisms are present at the same time and try to quantify the importance of these in natural communities [@Leibold2006; @Martorell2014; @Vergnon2009; @Kalyuzhny2014]. The main point is to understand which species level traits are important for community dynamics and which ones are unimportant [@Matthews2014], and this is related to the scale of observation. The problems of pattern and scale are critical in ecology [@Chave2013; @Levin1992], because patterns that seem stochastic at one scale may reveal structure at another scale. The concept of pattern is related to some sort of repetition that our brain can detect, but when this pattern repeats at different scales we talk about scale invariance or auto-similarity, and these patterns could be produced by critical phase transitions. Critical phase transitions were first introduced in ecology in the framework of habitat fragmentation [@Bascompte1996]. This phenomenon is characterized by the presence of two phases defined by some macroscopic features, that are linked by a critical point were a sudden transition happens [@Sole2006]. 

Several different ecological spatial models exhibit critical behavior related to the degree of disturbance [@Pascual2005a]. Some of these models showed robust criticality: a particular kind of criticality discovered for ecological systems [@Roy2003], where the scaling laws are present for a wide range of parameters. More important is that this kind of criticality has been documented for arid ecosystems [@Sole2007]; here the sudden shift towards a desert condition might occur when rainfall decrease [@Scanlon2007]or also with more intense grazing [@Kefi2007b]. The main mechanism is the positive effect produced by local facilitation, the chance of a new seedling to become established is higher near the parent plant.

Another example of an ecosystem exhibiting criticality are savannas, where the transition occurs between tree and grass cover [@Abades2014]. In critical phenomena the transition is produced by the capacity of the system to transmit some signal or information, in savannas the proportion of 60% grass 40% trees is linked to the threshold needed for fire to spread. The increase in the proportion of trees, due to a change in environmental conditions, can create positive feedback mechanism resulting in the encroachment of savanna ecosystems [@Abades2014].  

* NEUTRAL MODELS produce power laws [@Manor2008]!  

* Phase transition not linked to disturbances but to competence!!

(
Kalyuzhny M, Seri E, Chocron R, Flather CH, Kadmon R, Shnerb NM
Niche versus neutrality: a dynamical analysis.
The American naturalist. 

2014Fort H, Inchausti P Tropical forests are non-equilibrium ecosystems governed by interspecific competition based on universal 1/6 niche width.
PloS one. 2013
Condit R, Chisholm RA, Hubbell SP Thirty years of forest census at Barro Colorado and the importance of immigration in maintaining diversity.
PloS one. 2012 

)

Here I study the simplest model of neutral-niche communities where niche dynamics is represented as a competitive hierarchy. The competition intensity is represented as a parameter that modulates the transition between these extremes. 

Following this view I present a model that unify the Tilman's model of hierarchical competition with the classical neutral model using one parameter $\rho$ to link them in a spatially extended system. In this model a phase transition that produce a biodiversity collapse is observed. 

I use a stochastic cellular automata to show that there is a phase transition between the neutral and niche model with a spanning cluster formed by the most abundant species.

* biodiversity collapse [@Sole2004c]

* Multispecies early indicators [@Dakos2014]


* Communities evolve to avoid competition (Ghost of competition) thus they naturally have a very low competitive intensity. The effect of the niche is to keep low the competitive intensity 

* Fragmentation push the critical value to a lower competitive intensity - Comparison of fragments of different sizes with neutral dynamics rho=0 and with rho near critical value. Fragmentation makes the community to cross the critical value and become more hierarchical. 

See [@Pueyo2007] 

dos Santos FAS, Johst K, Grimm V (2011) Neutral communities may lead to decreasing diversity-disturbance relationships: insights from a generic simulation model. Ecol Lett 14: 653–660. doi:10.1111/j.1461-0248.2011.01626.x.

Reconciling neutral community models and environmental filtering: theory and an empirical test [@Jabot2008; ]

Quantifying the importance of local niche-based and stochastic processes to tropical tree community assembly [@Shipley2011]

When can we distinguish between neutral and non-neutral processes in community dynamics under ecological drift? [@Ruokolainen2009] 

Beta-Diversity in Tropical Forest Trees [@Condit2002; @Chave2002a]  

# Methods

## The spatial stochastic model

This model represent a continuum between hierarchical and neutral model in the same spirit as in [@Gravel2006;@Chisholm2010;@Zhou2008]. The model is a stochastic cellular automata (CA) or also called interactive particle system [@Durrett1994a]. In these kind of models space is discretized into a grid and only one individual can occupy a particular position. Each position represents an area fixed by the investigator to mimic the real system. Time is continuous so the update of the model is asynchronous. I update one randomly chosen site at a time and to perform one complete time interval $c J$ sites have to be updated, where $c$ is a constant that describes the overall rate at which transitions are occurring and $J$ is the size of the grid [@Durrett1994a].   

The model use periodic boundary conditions, which makes the landscape a torus. It means that sites on the top edge of the grid are neighbors of those on the bottom edge, and sites on the right edge are neighbors of those on the left. With this choice I can avoid edge effects and is equivalent to thinking that the grid is embedded in a large community.

The size of the community is given by *J = dimX* x *dimY*, where *dimX* and *dimY* are the dimension of the grid. Thus *J* is the maximum number of individuals in the simulated area. 

In this model all individuals have the same parameters, besides they should belong to different species [@Hubbell2001], and each species is assigned with a number. There are only two possible differences between species: 

* They may have a different frequency in the metacommunity and also different abundances in the local community.

* Hierarchical competition: species with lower numbers have a probability to replace species with higher numbers as in [@Tilman1994]. Thus a species with number 1 have a probability to replace species with number 2 and greater. The species with number 2 can replace species starting from 3. The probability of replacement is a parameter, when it is 0 replacement occurs only when a species dies.  

The colonization-competition and other trade-off are not explicitly included in the model. But a colonization-competition trade-off can be established if species numbering is arranged in inverse order as it's abundance $X_i$ in the metacommunity, the most competitive species (with number 1) will have the lowest migration rate and the less competitive will have the highest migration rate.  

There are four processes included in the model: death, local dispersal, and migration, starting with an empty site the following events can happen:

(1) With probability *m* an individual of a species *i* can migrate from the metacommunity at a rate proportional to its frequency $X_i$ in the metacommunity.

(2) When the grid is not full, individuals give birth with rate 1 to a new individual that disperse to the neighborhood with a dispersal kernel, here I use an inverse power kernel [@Marco2011].

(2) Individuals die a rate $\mu$

(5) When an individual dies it is replaced by a migrant from metacommunity with probability $m$ and with probability $1-m$ by an individual from the neighborhood. The neighborhood is established using the dispersal kernel with average distance $d$. Once the grid is full it keeps full, because when an individual dies is immediately replaced by another. This is called the zero-sum assumption in neutral models. 

(6) If the individual does not die it can be replaced by an individual from the metacommunity or neighborhood as in (5), but an individual of species with number $k$ can replace and individual of a species $k+1$ with probability $\rho$. Thus a hierarchical ordering of species is established. When this probability is zero the model behavior is neutral.

## Simulations

* I defined as tuning parameter the probability of replacement of individuals of different species $\rho$, the order parameter is  spanning cluster probability $SC_p$.

* The size of the lattice affects the value critical probability $p_c$ at which the transition
occurs; in small lattices $SC_p$ is non-zero for values of $\rho$ below the $p_c$, thus it is easier to form patches that connect the entire lattice. Therefore, in order to obtain an asymptotic estimate for the $p_c$ I performed a finite size scaling analysis. For this, I run simulations for different lattice sizes (Side = 100, 256, 512) and obtained asymptotic $p_c$ values by regressing $p_c$ against
$1/N$, the intercept becomes an estimate for a lattice of infinite size [@Stauffer1994; @Sornette2013].

* I determined critical probabilities for two different metacommunities: a) One with a logseries species abundance distribution, the most common distribution that fits experimental data [@White2012]. b) A uniform species distribution, this is analogous to simulate the apparition of a new species by evolution.    


# Results


# Discussion

* This phase transition is characteristic of some ecological systems and has been termed robust criticality. Here we describe for the first time the observation of this kind of criticality for a community assembly model. 


# References
