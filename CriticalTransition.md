# A spatial phase transition between neutral and niche ecological communities

**Leonardo A. Saravia** ^1^, **Jordi Bascompte** ^2^

1. Instituto de Ciencias Básicas

	Universidad Nacional de General Sarmiento

	J.M. Gutierrez 1159 (1613), Los Polvorines

	Buenos Aires, Argentina.

	<lsaravia@ungs.edu.ar>

2. Institute of Evolutionary Biology and Environmental Studies University of Zurich
	
	Winterthurerstrasse 190, 8057 Zurich Switzerland

	<jordi.bascompte@ieu.uzh.ch>

## Abstract

The dynamics of ecological communities can be described by two contrasting models: the first assumes that the individuals of all species are identical and do not have competitive interactions, the second assumes that species are different, adapted to particular habitat conditions and have strong interactions. These represent extremes of a continuum: the first is the neutral and the later the niche model of communities. Real communities are actually a mixture of both dynamics. Here I study the simplest model of neutral-niche communities where niche dynamics is represented as a competitive hierarchy. The competition intensity is represented as a parameter that modulates the transition between these extremes. I use a stochastic cellular automata to show that there is a phase transition between the neutral and niche model with a spanning cluster formed by the most abundant species. The transition implies a sharp fall of species diversity (Shannon entropy) but the richness (number of species) shows a gentle decline with increasing competitive intensity. The critical point is at a very low value of competitive intensity and the same was reported for different real communities. Thus it will be possible that a mechanism of self-organization is driving communities to criticality, that mechanism could be similar but not equal to self-organizing criticality.


## Introduction

Much effort has been devoted to understand the mechanisms of community assembly and dynamics. In principle the emphasis were on deterministic mechanisms based on niche differences between species; the niche theory assumes that different species are regulated by different environmental factors and infer that diversity originates from spatial and temporal environmental heterogeneity [@Hutchinson1957; @Tilman1982; @Chesson2000]. More recently the emphasis shifted to stochastic mechanisms in the form of the Neutral theory of Biodiversity and Biogeography [@Hubbell2001]. The neutral theory assumes that individuals of all species are equivalent and it proposes that diversity originates from a balance between immigration, speciation, and extinction. Neutral theory has been proposed as a parsimonious formulation that can provide new insight into the patterns of community assembly [@Hubbell2005], besides this simplification it can predict some community metrics very well [@Volkov2007; @Rosindell2012], mainly the species abundance distribution (SAD).

Finally arises a unified view that accepts that both kinds of mechanisms are present at the same time and try to quantify the importance of these in natural communities [@Leibold2006; @Martorell2014; @Vergnon2009; @Kalyuzhny2014]. The main point is to understand which species level traits are important for community dynamics and which ones are unimportant [@Matthews2014], and this is related to the scale of observation. The problems of pattern and scale are critical in ecology [@Chave2013; @Levin1992], because patterns that seem stochastic at one scale may reveal structure at another scale. The concept of pattern is related to some sort of repetition that our brain can detect, but when this pattern repeats at different scales we talk about scale invariance or auto-similarity, described by power laws. These patterns could be produced by critical phase transitions described by percolation theory [@Stauffer1994] . Critical phase transitions were first introduced in ecology in the framework of landscape ecology [@Loehle1996b] and habitat fragmentation [@Bascompte1996]. Percolation is characterized by the presence of two phases defined by some macroscopic features, that are linked by a critical point were a sudden transition happens  and a large spatial pattern emerges. In a two dimensional landscape where each site is connected to the others neighbors with some probability $p$, if you start with a small $p$ there will be some connected clusters, but if you increase $p$ there is a point where a cluster of sites spanning the entire landscape appears. This spanning cluster has a self-similar structure and is produced by local interactions [@Sole2006]. 

Several different ecological spatial models exhibit critical behavior related to the degree of disturbance [@Pascual2005a]. Some of these models showed robust criticality: a particular kind of criticality discovered for ecological systems [@Roy2003], where the scaling laws are present for a wide range of parameters. More important is that this kind of criticality has been documented for arid ecosystems [@Sole2007]; here the sudden shift towards a desert condition might occur when rainfall decrease [@Scanlon2007] or also with more intense grazing [@Kefi2007b]. The main mechanism is the positive effect produced by local facilitation, the chance of a new seedling to become established is higher near the parent plant.

(Describe SOC as another possible mechanism?)

Another example of an ecosystem exhibiting criticality are savannas, where the transition occurs between tree and grass cover [@Abades2014]. In critical phenomena the transition is produced by the capacity of the system to transmit some signal or information, in savannas the proportion of 60% grass 40% trees is linked to the threshold needed for fire to spread. The increase in the proportion of trees, due to a change in environmental conditions, can create positive feedback mechanisms resulting in the encroachment of savanna ecosystems [@Abades2014].  

Neutral models can produce species patches with power law distributions without been near a critical state [@Houchmandzadeh2003], so the detection of this kind of patterns do not imply a phase transition. Moreover phase transitions are also observed in non-spatial models, some studies detected the presence of a sharp transition between neutral and niche dynamics in spatially implicit models [@Chisholm2011] and Fisher et al. [-@Fisher2014] demonstrated the presence of a phase transition for this kind of neutral-niche models.

Here I study a different kind of ecological spatial phase transition that is not related to disturbance or fire, is the transition between a neutral and a niche community. I formulated the simplest model of neutral-niche communities where niche dynamics is represented as a competitive hierarchy [@Saravia2015]. This spatially explicit model unifies the Tilman's model of hierarchical competition with the classical neutral model using one parameter: the competition intensity. This parameter is represented as the probability that one species replace another and modulates the transition between the neutral phase and niche phase. 

Our first aim is to demonstrate the existence of a spatial phase transition in neutral-niche models, our second objective is to explore the dependence of the critical point with the parameters that are more important for the spatial pattern and to suggest some possible early warnings of the transition. 


# Methods

First we define the spatial explicit neutral-hierarchical model, then we explain how we characterized its critical behavior in terms of percolation theory and how we made simulations. We refer to interested readers to more extensive introductions to percolation theory in an ecological context [@Sole2006;Oborny2007]. 

## The spatial stochastic model

This model represent a continuum between hierarchical and neutral model in the same spirit as in [@Gravel2006;@Chisholm2010;@Zhou2008]. The model is a stochastic cellular automata (CA) or also called interactive particle system [@Durrett1994a]. In these kind of models space is discretized into a grid and only one individual can occupy a particular position. Each position represents an area fixed by the investigator to mimic the real system. Time is continuous so the update of the model is asynchronous. I update one randomly chosen site at a time and to perform one complete time interval $c J$ sites have to be updated, where $c$ is a constant that describes the overall rate at which transitions are occurring and $J$ is the size of the grid [@Durrett1994a].   

We use periodic boundary conditions, which makes the landscape a torus. It means that sites on the top edge of the grid are neighbors of those on the bottom edge, and sites on the right edge are neighbors of those on the left. With this choice we can avoid edge effects and is equivalent to thinking that the grid is embedded in a large community.

The size of the community is given by *J = dimX* x *dimY*, where *dimX* and *dimY* are the dimension of the grid. Thus *J* is the maximum number of individuals in the simulated area. As in a classical neutral model there is a metacommunity i.e. a regional species pool assumed very large and invariant in ecological time scales [@Hubbell2001]. All individuals have the same parameters, besides they should belong to different species [@Hubbell2001], and each species is assigned with a number. There are only two possible differences between species: 

* They may have a different frequency $X_i$ in the metacommunity and also different abundances in the local community.

* Hierarchical competition: species with lower numbers have a probability to replace species with higher numbers as in [@Tilman1994]. Thus a species with number 1 have a probability to replace species with number 2 and greater. The species with number 2 can replace species starting from 3. The probability of replacement is a parameter, when it is 0 replacement occurs only when a species dies.  

The colonization-competition and other possible trade-offs are not explicitly included in the model. But a colonization-competition trade-off can be established if species numbering is arranged in inverse order as it's frequency $X_i$ in the metacommunity, the most competitive species (with number 1) will have the lowest migration rate and the less competitive will have the highest migration rate.  

There are four processes included in the model: death, local dispersal, and migration, starting with an empty site the following events can happen:

(1) With probability *m* an individual of a species *i* can migrate from the metacommunity at a rate proportional to its frequency $X_i$ in the metacommunity.

(2) When the grid is not full, individuals give birth with rate 1 to a new individual that disperse to the neighborhood with a dispersal kernel, here I use an inverse power kernel [@Marco2011]:

 	$d(x) =  \frac{\alpha -1}{x_{min}} \left(\frac{x}{x_{min}} \right)^{-\alpha}$ with $mean =\frac{\alpha-1}{\alpha-2}x_{min}$ where $\alpha > 1$ and $x \ge x_{min}$. 

	where $d(x)$ is the probability that an individual disperse a distance $x$ from the parent. In all cases I used $x_{min} = 1$.


(2) Individuals die a rate $\mu$

(5) When an individual dies it is replaced by a migrant from metacommunity with probability $m$ and with probability $1-m$ by an individual from the neighborhood. The neighborhood is established using the dispersal kernel with average distance $d$. Once the grid is full it keeps full, because when an individual dies is immediately replaced by another. This is called the zero-sum assumption in neutral models. 

(6) If the individual does not die it can be replaced by an individual from the metacommunity or neighborhood as in (5), but an individual of species with number $k$ can replace and individual of a species $k+1$ with probability $\rho$. Thus a hierarchical ordering of species is established. When this probability is zero the model behavior is neutral.

The model was developed using the C++ programing language and its source code is available at <https://github/lasaravia/neutral> and figshare <http://dx.doi.org/10.6084/m9.figshare.969692>. 

## Percolation and simulations

To characterize our model in terms of percolation theory we need to define an order parameter that depends on an external control or tunning parameter that can be continuously varied. Thus we defined as tunning parameter the replacement probability $\rho$, and the order parameter was the spanning cluster probability $SC_p$. 

In our model percolation is produced when there is at least one patch of one species that spans from one edge of the system to the opposite edge. To detect species patches we used a modified Hoshen–Kopelman cluster labeling algorithm [@Hoshen1976] with a neighborhood defined by the four nearest sites (Von Neumann neighborhood) available at github (https://github.com/lsaravia/Clusters). The percolation point is defined as the  value of the tunning parameter $\rho$ at which $SC_p$ is 0.5 and the patch size distributions were measured after 5000 time iterations. 

The size of the lattice affects the value of the critical probability $p_c$ at which the transition
occurs; in small lattices $SC_p$ is non-zero for values of $\rho$ below the $p_c$, thus it is easier to form patches that connect the entire lattice. Therefore, in order to obtain an asymptotic estimate for the $p_c$ we performed a finite size scaling analysis. For this, we run simulations for different lattice sizes (Side = 100, 150, 256, 512) and obtained asymptotic $p_c$ values by regressing $p_c$ against
$1/N$, the intercept becomes an estimate for a lattice of infinite size [@Stauffer1994; @Sornette2013].

We determined critical probabilities for two different metacommunities: a) One with a logseries species abundance distribution, the most common distribution that fits experimental data [@White2012]. With this metacommunity we included a competition-colonization trade-off by arranging species numbers in inverse order as it's frequency $X_i$ in the metacommunity. b) A uniform species distribution, all species have the same probability to colonize the local community, this is analogous to simulate the apparition of a new species by evolution thus the migration parameter $m$ would represent an speciation rate. The values of the $m$ parameter (Table 1) were at least two orders of magnitude higher than realist speciation rates [@Rosindell2009] thus we could not interpret our results in an evolutionary framework.    

All simulations started with an empty lattice that is colonized by migrants mimicking the assembly of a new community. The range of parameters used were compatible with published results in tropical forest [@Anand2010; @Condit2002; @Etienne2007], It was suggested that fat-tail dispersal kernels give more realistic results [@Rosindell2009; @Seri2012] so we used an inverse power law distribution with an exponents always greater than two---so the mean exist (Table 1). The parameter $\rho$ is varied in all the range between 0 and 1 to determine the critical point, in the region where we suspect the $p_c$ will be located (near 0) the steps were very small (0.0001) and greater (0.1) in the region were we don't expect the $p_c$. All the simulations and the analysis of the model output were done in the R statistical statistical language [@RCoreTeam2014] and the scripts are available at github <https://github.com/lsaravia/CriticalTransition>




# Results


# Discussion

* This phase transition is characteristic of some ecological systems and has been termed robust criticality. Here we describe for the first time the observation of this kind of criticality for a community assembly model. 

* Communities evolve to avoid competition (Ghost of competition) thus they naturally have a very low competitive intensity. The effect of the niche is to keep low the competitive intensity. This suggest the possible existence of mechanism similar to critical self organization. 

* Fragmentation push the critical value to a lower competitive intensity - Comparison of fragments of different sizes with neutral dynamics $\rho=0$ and with rho near critical value. Fragmentation makes the community to cross the critical value and become more hierarchical. 

* The existence of an spanning patch of one dominant species can point us that the system is near the critical point.

* This implies that as spatially implicit models suggest [@Chisholm2011] most ecosystems will exhibit patterns of diversity that are either strongly niche-structured or indistinguishable from neutral.

* Fisher2014 after disturbances their model predicts that under stress a community will suffer a biodiversity collapse produced by a shift towards neutrality, this means that disturbed less diverse communities should have neutral dynamics. This prediction is contrary to most of the models and experimental data that suggest that niche dynamics dominate low-diversity communities while neutral dynamics will be more common in high diversity communities  [@Chisholm2011]


* At the smaller lattice sizes (sides 100-150 in our case) there is a non-zero probability of a spanning clusters ($SC_p$) in the neutral phase, but at after that the $SC_p$ drops to zero. Abades et al. [-@Abades2014] reported a critical phase transition for a neutral model using a lattice with sides from 10 to 100, their results are analogous to our neutral phase. Only with lattice sides lower than 150 we observed a non-zero probability of a spanning clusters ($SC_p$), after that the $SC_p$ drops to zero. Besides that they use a different tunning parameter and a different neighborhood there is some concerns that their results could be due to the limited range of size they used in the spatial simulations. 

* Detections of communities at critical state:
 	* long transient time 
 	* spanning patch one species occupying 30% of the area..

* New directions
	* Study realistic evolutionary rates m
	* Study how the diversity of metacommunity influence
	* Disturbances!
	* Extension to food webs trophics


# References
