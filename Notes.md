# Notes

### Ecology reviewer 1

More model repetitions to estimate critical point--> 100?

More repetitions to study de patch distributions --> 100?

5000 iterations are enough to get a steady state? Add an algorithm to estimate steady states using H?

Using exponential decay for early warning transitions does not convice (Also JBonachela)  

The exponential tail being may not be a result of the finite space in the simulation, so better not mention it as an e

Delete quantil regressions

### Para agregar a la discusion

Ni los tradeoff ni el dispersal pueden cambiar la sensibilidad a la intensidad de la competencia, por lo tanto lo mas probable es que la intensidad de la competencia en comunidades naturales sea muy baja.


### percolation

The simplest mechanism that produces scale invariance is percolation, physically describes the movement of a fluid through porous media. One can think about a two dimensional landscape where each site is connected to the others neighbors with some probability $p$. The landscape will be connected when $p$ reaches some threshold value . 


The most common species in a neutral community should coincide with the most common species in the metacommunity, probably a species with high dispersal capabilities. But a hierarchycal community, if the community is at steady state, the most common species will be the most competitive. If you study the traits of the most common and the community is near steady state you will know what are the principal forces structuring the community.


* biodiversity collapse [@Sole2004c]

* Multispecies early indicators [@Dakos2014]

* The most competitive species could be seem as non-habitable patch produced by deforestation, and the intensity of competition is a measure of growth rate of deforestation. It has been seen that deforested patches become abandoned due to economic crisis or changes reglamentations and management [@Haddad2015;Hansen2013]. Besides the growth and decay of unhabitable patches should be different from other species rates of invasions (here we assumed all equal) the qualitative features should be the same taking into account that near the percolation threshold the details of the models do not modify the behavior of the system. Thus based on these we can establish that about a 30% percent non-habitable patches can produce a biodiversity collapse if communities live near a critical interaction/competitive intensity. A similar approach taking into account the interplay between two dynamic processes: habitat loss/recuperation and 1 species population spread [@Oborny2007]. found a critical threshold on the same order: 27% of non-habitable patches. 

* Evolution: we could start with a neutral community and let species mutate to modify its interaction by individual ver Eco-evolutionary dynamics

## Reviewers

Fangliang He,

Stephen P. Hubbell 
shubbell@eeb.ucla.edu
UCLA Ecol & Evolutionary Bio

Samir Suweis
samir.suweis@pd.infn.it
 Dipartimento di Fisica e Astronomia

Carlos Martorell 
martorell@ciencias.unam.mx ,
Facultad de Ciencias, Departamento de Ecología y Recursos Naturales, Universidad Nacional Aut o
Circuito Exterior S/N, Cd. Universitaria, 04510 México D.F., Mexico; 

James Rosindell 
j.rosindell@imperial.ac.uk
Faculty of Natural Sciences, Department of Life Sciences (Silwood Park)


C. Patrick Doncaster
cpd@soton.ac.uk
School of Biological Sciences, University of Southampton, Southampton, United Kingdom


# Fragments from papers



-------------------

1. Ostling AM (2012) Large-scale spatial synchrony and the stability of forest biodiversity revisited. Journal of Plant Ecology 5: 52–63. doi:10.1093/jpe/rtr035.

I used a mean dispersal distance r = 40 m, which is the average r across species of fits to seed trap data for a tropical forest plot in Panama (Condit et al. 2002). Based

------------------

1. Condit R, Pitman N, Jr. EGL, Chave J, Terborgh J, et al. (2002) Beta-Diversity in Tropical Forest Trees. Science (80- ) 295: 666–669.

Beta-diversity:

The Sørensen index is $S_12/[0.5(S_1 + S_2 )]$, where S12
is the number of species common to both sites and
Si is the total found at site i. Jaccard’s index is
S12/(S1 + S2 - S12 ). For two plots, the probability
F can be calculated as F = SUM ͚f_i1 f_i2, where fij is the
relative abundance of species i at site j, and the
sum is over all species at both sites. More generally, F(r) can be calculated by finding all pairs of
trees separated by distances between r and $r + \Delta r$
and then determining what proportion of these
pairs are the same species.


-----------------


May RM, Crawley MJ, Sugihara G (2007) Communities: patterns. 
In: May RM, McLean AR, editors. Theoretical ecology: principles 
and applications. New York: Oxford University Press. pp. 111–131. 

pag 119 

A serious difficulty with this neutral model is its extreme sensitivity 
to the neutrality assumption. 

Once differences between species are allowed, then competitive 
exclusion is likely to occur rapidly. 

-----------

More generally, the June 2006 issue of Ecology contained a series 
of papers in a special feature on 

neutral community ecology (Naeem, 2006). 

The degree of apparent community structure—niche rather than 
neutral— depends on many factors, including number of species (diversity), 
degree of niche overlap, dispersal abilities, and assumptions 
made about environmental variability (Gravel et al., 2006). 

The theme is re-echoed, with variations, by Bonsall et al. (2004b), 
who show how ‘organization dominated by niche structure and organization through chance and 
neutral processes’ can operate simultaneously, with varying weights, depending on the interplay between population dynamics (ecology) and life-history trade offs (evolution). 

(Ver Bonsall, M.B., Jansen, V.A.A., and Hassell, M.P. 2004b. 

Life history trade-offs assemble ecological guilds. Science 
306: 111–114). 

page 131 

As in all such quests, it is good to keep in mind Einstein’s dictum, 
which translates roughly as ‘seek simplicity, and distrust 
it’. 



-----------

Chisholm RA, Pacala SW (2010) Niche and neutral models predict 
asymptotically equivalent species abundance distributions 
in high-diversity ecological communities. PNAS 107: 15821–15825. 
doi:10.1073/pnas.1009387107. 

Some ecologists feel uncomfortable with neutral theory be cause 
it is wrong in the sense that it assumes species equivalence in 
the face of a vast literature cataloguing species differences 
and because it subsumes many undoubtedly real and well-researched 
ecological processes under the umbrella of stochasticity, 
which does not itself really constitute a force of nature (26). 

----------------

Dakos2014

We find that critical slowing-down indicators derived from time series of biomasses measured at the species and community level signal the proximity to the onset of community collapse. In particular, we identify specialist species as likely the best-indicator species for monitoring the proximity of a community to collapse. In addition, trends in slowing-down indicators are strongly correlated to the timing of species extinction.

----------------


@Keymer1998

Fixed-point attractor: we say that the dynamic
converged to a fixed-point attractor if the system
did not change in the frequency of occupied
patches in two successive generations.
Statistical steady-state: we say that the
dynamic converged to a statistical steady-state
when the proportion of cells in each state, as
measured by the Shannon entropy of the system,
remain bounded to very small fluctuations
around a fixed value, although the system
changes with no apparent trend over time.
Non-convergent: we say that the system did
not converge if after 10 500 generations we do
not find any of the above mentioned attractors.

---------

@Gravel2006 Reconciling niche and neutrality: the continuum hypothesis.


For simplicity, we assume species are equal in their death probability, although it is obvious that differentiation could also occur through this process. The models also focus on the local community scale. Thus, we do not consider speciation, and we recognize that further work is necessary to unify the models at the metacommunity scale.

	
It has been argued that observing patterns of species abundance 
typically predicted by neutral models cannot necessarily be 
taken as evidence of neutral processes (Purves & Pacala 2005), 
and our results concur with that. However, our results also show 
that the observed patterns can result from a neutral drift created 
by an elevated niche overlap, sustained by the constant reintro- 
duction of excluded species through high immigration and/ or 
speciation rates. Under such conditions, the rare species are 
transients and their occurrence depends on immigration; on 
the other hand, the resident species are permanent and their 
occurrence depends on niche differentiation (Magurran & Henderson 
2003; Schwilk & Ackerly 2005). 


Results demonstrate that niche and neutrality form ends of a continuum from competitive to stochastic exclusion. In the absence of immigration, competitive exclusion tends to create a regular spacing of niches. However, immigration prevents the establishment of a limiting similarity.

A community solely driven by competition will have a
deterministic succession, while a neutral community will
have the maximum stochasticity among replicated succes-
sions because of the neutral drift (Clark & McLachlan 2003).
We used this proxy to assess ÔneutralityÕ when we tested the
different predictions presented above.

The index was calculated after 500 time steps. This is not a sufficient
run length to reach a truly stable species composition, but
preliminary simulations showed that index differences
among scenarios were similar when the run length was
longer than 500 time steps

-----------------------
@Kefi2011

The percolation probability was
estimated using 500 · 500 lattices. The system is considered to be percolated when it has
at least one patch that spans from one edge of the system to the opposite edge. The
location of the percolation point is defined as the parameter value at which the percolation
probability is 0.5 (see the Results for further details about the percolation point).
The patch size distributions were obtained as follows. Simulations were run for
10 000 time units. The first 3000 transient time steps were discarded. Thereafter, patch
size data were recorded every 40 snapshots, to minimize any temporal correlation
among successive snapshots. Data from 175 such snapshots were used to plot the patch
size distribution.

# discussion
Our results suggest that tracking changes in patch size distributions
associated with knowledge on the connectedness of the system (and more specifically
on the presence or absence of spanning clusters) might provide qualitative information
regarding the level of stress exerted on the ecosystem.

Additionally, the patch size distribution can provide some information about
underlying ecological mechanisms. We showed that systems with the same cover (but
different underlying ecological mechanism) can be described by different distributions
(Fig. 2; Appendix S3). Considering a real ecosystem and knowing its fractional cover,
the patch size distribution can be compared with the patch size distribution generated
by a null model for the same fractional cover. An excess in large patches compared with
the null model would mean that local interactions might be at play, as a mechanism
promoting longer spatial correlations than expected otherwise.

-----------------------

@Orrock2010 Local community size mediates ecological drift and competition in metacommunities

The model highlights the role of community size, i.e. the number of competitors in the local community, in mediating the relative importance of stochastic and deterministic forces. In metacommunities where local communities are small, ecological drift is substantial enough that strong competitors become effectively neutral, creating abrupt changes in the outcome of competition not predicted by the standard competition-colonization trade-off. Importantly, the model illustrates that, even when other aspects of species interactions (e.g. migration ability, competitive ability) are unchanged, local community size can alter the dynamics of metacommunity persistence. Our work demonstrates that activities which reduce the size of local communities, such as habitat destruction and degradation, effectively compound the extinction debt

------------------------
@Fort2013a  Tropical Forests Are Non-Equilibrium Ecosystems Governed by Interspecific Competition Based on Universal 1/6 Niche Width

This niche width yields an average ratio of 0.25 between interspecific and intraspecific competition that corresponds to an intermediate value between the extreme claims of the neutral model and the classical niche-based model of community assembly (where interspecific competition is dominant).


------------------------
@Martorell2014 Testing the roles of competition, facilitation and stochasticity on community structure in a species-rich assemblage

Most ecologists would agree that interactions and stochasticity occur concurrently in ecological communities. Neutrality may operate within subsets of closely related or convergent species (Leibold & McPeek 2006; Scheffer & van Nes 2006; Vergnon, Dulvy & Freckleton 2009), or among species pairs with similar fitness (equalizing mechanisms). However, increasing niche differentiation is required for coexistence to occur as asymmetries in competitive ability become greater (stabilizing mechanisms; Chesson 2000b; Adler, HilleRis- Lambers & Levine 2007). Furthermore, immigration may counteract competitive exclusion in interconnected metacommunities, reducing the degree of niche organization required for diversity maintenance (Gravel et al. 2006). Nevertheless, stochasticity is expected to set more restrictive conditions to coexistence in interaction-driven communities (May 1973; Tilman 2004; Adler & Drake 2008).

Empirical evidence is still insufficient to settle the conflict between these theories. Community patterns expected from theories based on interactions and stochasticity are frequently indistinguishable (Hubbell 2001; Chave, Muller-Landau & Levin 2002; Rosindell, Hubbell & Etienne 2011),

The fact that intraspecific competition (altogether with stochasticity) was able to explain
most of the variance in species abundance in our community would appear to support his claims. Nevertheless, the same pattern would be expected if the community were structured by interspecific competition, which would cause extinctions until only a set of species that are more limited by conspecific than heterospecific interactions remained (Chesson 2000b).

Our results, taken together with recent studies on tropical forests, suggest that weak interactions among established plants may be a general phenomenon, but that local interactions during colonization are important drivers of community composition.

----------------------
@Volkov2009

Both approaches yield very similar answers and show
that the collective effects of the pairwise interspecific interaction
strengths are weak compared with the intraspecific interactions.
Our approaches can be applied to other ecological communities
in steady state to evaluate the extent to which interactions need
to be incorporated into theoretical explanations for their structure
and dynamics.

However, both models have different underlying assumptions and
complement each other, and the degree of accord is noteworthy.
In general, from a quantitative perspective, both approaches lead
to the conclusion that the effects of interspecific interactions are
relatively weak and impose a relatively minor constraint on species
dynamics in the BCI forest (Fig. 3). It is as though evolution has
chosen weakly interacting species for proximal coexistence.


----------------------
@Vergnon2009 Niches versus neutrality: uncovering the drivers of diversity in a species-rich community

Our study strongly suggests that a combination of niche-based and equalizing processes best explains the patterns and species dynamics within the L4 station marine phytoplankton community (Table 2).


---------
@Xu2015a Partial recovery of a tropical rainforest a half century after clear-cut and
selective logging

Species richness recovered faster than species composition and structure in both selectively
logged and clear-cut forests. Both total number of species and number of rare spe

--------

@Chisholm2010 Niche and neutral models predict asymptotically equivalent species abundance distributions in high-diversity ecological communities

niche theory and neutral theory are two major families of theoretical models that aim to explain patterns of biodiversity observed in nature (1, 2). Niche theory, which has a long history of development (3–5), assumes that different species are regulated by different environmental factors and proposes that diversity arises from spatial and temporal environmental heterogeneity (6). Neutral theory, which was developed more recently from the theory of island biogeography (7), assumes that species are equivalent and proposes that diversity arises from a balance between immigration, speciation, and extinction (8–10)

in which a semiisolated local community receives immigrants from a much larger metacommunity that is assumed to be static relative to the timescale of the local community 

Our results confirm that the neutral model of species abundances emerges from more complex models (29) at sufficiently large scales, sufficiently high levels of diversity, and sufficiently low niche dimensionality.

-------
Eco-evolutionary responses of biodiversity to climate change.

F1000 Recommendation - DOI: 10.3410/f.717967699.793467371

This interesting study investigates the role of multi-species interactions and evolutionary processes in determining species extinction due to a changing climate. Current predictions of the effect of climate change on biodiversity commonly overlook processes of species interactions and the possibility of adaptation to a changing climate. In this study, Norberg et al. develop a spatially explicit eco-evolutionary model of multi-species responses to climate change. Their study suggests that species with higher genetic variance and shorter dispersal distances are less likely to become extinct due to climate change. This is because reduced connectivity facilitates the fixation of mutations that may allow species to adapt to the new environmental conditions. 

This contradicts the current understanding that populations well-connected at large spatial scales are less vulnerable to disturbances. The rationale behind this idea is that isolated populations have a greater difficulty to recover from disturbances, since they do not receive recruits from other locations. In terms of management for conservation, this has led to an effort in designing networks of protected areas that maximise the connectivity among populations (e.g. spacing of marine protected areas and green corridors). While maximising connectivity is a valid practice to minimise risks of extinction due to habitat fragmentation and episodic acute disturbances, Norberg et al.'s findings suggest it might not be ideal to face climate change. Thus, instead of aiming to 'maximise' connectivity, maybe we should start aiming to 'optimise' connectivity. 

The authors also highlight the role of competition in the extinction of species and predict the occurrence of extinction and evolutionary debts, i.e. changes in species richness that occur long after the climate stabilises.

----------
@Rosindell2011 The Unified Neutral Theory of Biodiversity and Biogeography at Age Ten.

Many different models can explain the same sampled
species abundance data [16–18]. Even models based on
neutral theory, but introducing a particular niche struc-
ture, give mathematically identical results to the classic
neutral model at large spatial scales [19,20]. On the one
hand, this implies that neutral models are robust to the
introduction of certain niche structures, but on the other
hand, it means that species abundance data will not always
be able to distinguish various models.


For example, fits of species–area curves by neutral models were not realistic [26,27] until long-distance dispersal was incorporated [28] (Box 3). This supports the general concept that long-distance dispersal is
important [29].

-------------
@Hubbell2005 Neutral theory in community ecology and the hypothesis of functional equivalence


The importance of trade-offs in resource use led Tilman (1994) to consider trade-offs in other life-history characteristics besides those involving resource exploitation, and here additional problems with the classical paradigm began to surface. In particular, Tilman revisited the old fugitive species concept of Hutchinson, and showed theoretically that if there was a strict transitive trade-off between competitive ability (site tenacity) and dispersal ability then, in principle, any arbitrary number of species could coexist.

However, soon thereafter, Hurtt & Pacala (1995) showed that Tilman’s trade-off assumption could be relaxed and that species of arbitrary competitive and dispersal abilities could coexist if dispersal and recruitment limitation were sufficiently strong. Hurtt and Pacala defined dispersal and recruitment limitation as the failure of species to reach and/or establish in all sites favourable for their growth and survival. When dispersal was infinite, such that all species reached all sites, then the best competitor won each site for which it was superior, the equilibrium solution expected from the classical theory. However, when dispersal and recruitment limitation were strong, many sites were won by default by competitively inferior species because the best competitor for the sites failed to reach them. Although the coexistence produced by this mechanism is non-equilibrium, Hurtt and Pacala showed that competitive exclusion could be essentially infinitely delayed. This result meant that Gause’s competitive exclusion principle was not generally true in spatial ecology. The slowness to reach equilibrium also depended on species richness.

---------------

@Jabot2008 Reconciling neutral community models and environmental filtering: theory and an empirical test

Value of m

--------------
@Condit2012 Thirty Years of Forest Census at Barro Colorado and the Importance of Immigration in Maintaining Diversity

Value of m

--------------
@Sole2007 Scaling laws in the drier

A fuller theory of scaling behaviour is required to provide firmer connections between predicted and observed patterns of desertification---and so to provide a better understanding of the nature of transitions between green and desert phases. These belong to the family of non-equilibrium phase transitions seen in several different fields of research8,9. There is a great opportunity here for interdisciplinary work that would have potentially far-reaching consequences in conservation biology.

Such power laws occur in other types of ecosystem and are a fingerprint of self-organization: that is, they are the result of internal dynamic processes driven by local interactions. This principle applies to the field data reported by both Scanlon et al.4 and Kéfi et al.5. It indicates that plant interactions play a central role in shaping these ecosystems, which as a whole are characterized by productivity levels that largely depend on precipitation.


---------------
@Abades2014 Fire, percolation thresholds and the savanna forest transition: a neutral model approach

Although many factors can promote increases in tree cover, acting alone or in tandem the 40% tree cover threshold is likely to be the tipping point for encroachment through disrupting grass percolation and thus suppressing fire spread

---------------
@Manor2008 Facilitation, competition, and vegetation patchiness: From scale free distribution to patterns

Recent reports (Scanlon et al., 2007; Kefi et al., 2007) have demonstrated power-law statistics of the patch size for vegetation ecosystems in the arid and semi-arid climatic zone across a wide range along the annual rainfall gradient. As pointed out in Scanlon et al. (2007) and Kefi et al. (2007), this phenomenon is actually puzzling.


Another mechanism that yields power-laws and other fat-
tailed distributions is multiplicative noise. This situation occurs
when the random fluctuations that affect the system are
proportional in magnitude to the size of the system itself. An
example is the neutral theory of species abundance (Hubbell,
2001); if the chances of any individual to produce an offspring and to die are the same, the abundance fluctuates along generations and the size of the fluctuation is proportional to the population.
This ‘‘law of proportion effect’’ was first discovered in the context
of business firm’s size (Gibrat, 1930; Simon and Bonini, 1958; Levy
and Solomon, 1996) and is also relevant to the effect of small Facilitation, competition, and vegetation patchiness: From scale free distribution to patterns
fluctuations in growing populations, like surname abundance
(Manrubia et al., 2003) and degree distribution in scale free networks (Barabasi et al., 1999).

---------------
@Stumpf2012 Critical Truths About Power Laws

Moreover, knowledge of whether or not a distribution is heavy-tailed is far more important than whether it can be fit using a power law.

Suppose that one generates a large number of independent random variables xi drawn from heavy-tailed distributions, which need not be power laws. Then, by a version of the central limit theorem (CLT), the sum of these random variables is generically power-law distributed (www.informs-sim.org/wsc04papers/016.pdf).



---------------

@Zillio2008 Incipient criticality in ecological communities

Our results demonstrate that scaling provides a model-independent
framework for analyzing and unifying ecological data and that,
despite the absence of power laws, ecosystems are poised in the
vicinity of a critical point.

----------------

@Fisher2014 The transition between the niche and neutral regimes in ecology


As a result, it is possible to drive a community between the niche and neutral phases by changing the environmental conditions. As an example, we consider the effects of selective logging on a population of butterflies in a tropical forest on Buru, Indonesia (32). Through habitat destruction, logging essentially moves the butterfly community from a position with high K/ω to one with low K/ω, tracing a path along the stochasticity axis in the phase diagram (Fig. 4A). Our model predicts that when a diverse community within the niche phase is placed under a stress that lowers K/ω to the critical value, it will undergo a transition to the neutral regime. LV sim- ulations show that this transition results in a collapse of biodiversity and leads to an increase in the skewness of the species abundance distribution 



-----------

@Rosindell2009 Species-area curves, neutral models, and long-distance dispersal

In the present manuscript, we extend to the case of ‘‘fat-tailed’’ dispersal, where the probability of dispersal decays as a power law with distance, which is considered a more faithful description for many species, especially wind-dispersed seeds (Clark et al. 1999, Nathan 2006)

Condit et al. (2002) quote a range of m values between 4.8 3 10À8 and 1.7 3 10À14 for forests in
Panama, Ecuador, and Peru,


Speciation rates are estimated at mean=0.56 SE .17 species
per million years (Baldwin and Sanderson 1998), which
corresponded to rates from 2.3 10^-6  to 1.1 10^-5 per
species per generation; note, however, that this figure
needs to be divided by the number of individuals per
species to give the speciation rate per birth m.




----------

@Seri2012 Neutral Dynamics and Cluster Statistics in a Tropical Forest.

The Cauchy kernel, with the right dispersal parameter
g (see “Methods”), was found to yield the observed pat-
terns with very good accuracy, while the data rule out the
entire spectrum of kernels covered by the MLGK. First,
we discuss the inadequacy of the MLGK (and the discrim-
inatory power of our technique), and then we show how
the Cauchy kernel yields a good fit to the data.

---------------------------
@Chase2007 Drought mediates the importance of stochastic community assembly

Here, I develop the hypothesis that the relative importance
of stochastic ecological drift and/or priority effects depend on the
harshness of the ecological filter in those habitats. I established
long-term experimental ponds to explore the relative importance
of community assembly history and drought on patterns of com-
munity compositional similarity among ponds that were otherwise
similar in their environmental conditions. I show considerable
site-to-site variation in pond community composition in the ab-
sence of drought that likely resulted from a combination of
stochastic ecological drift and priority effects. However, in ponds
that experienced drought, I found much higher similarity among
communities that likely resulted from niche-selection filtering out
species from the regional pool that could not tolerate such envi-
ronmental harshness


------------------------
@Kefi2014

If the patterns are non-periodic or irregular
(in particular, if the patch-size distribution is described by a heavy-
tailed distribution), the patch size distribution may contain
information about the degradation level of the ecosystem but
more needs to be known about the underlying ecological
mechanisms to interpret the changes in the shape of this
distribution [15,20,46,60].


There are other
promising avenues of research to be pursued. Composite metrics
that combine spatial patterns with their temporal dynamics could
potentially offer new and potentially more reliable indicators of
imminent transitions [66]

------------------------
@Newman2005

What does it mean to say that a distribution has no finite mean? Surely we can take the data for real solar flares and calculate their average? Indeed we can, but this is only because the data set is of finite size. When α ≤ 2, the value of the mean x is dominated by the largest of the samples drawn from the distribution. For any data set of finite size this largest sample has a finite value and hence so does the mean. But the more samples we draw, the larger, on average, will be the largest of them, and hence the larger the mean as well. The divergence of Eq. (11) is telling us that as we go to larger and larger data sets, our estimate of the mean x will increase without bound. We discuss this more below.


-------------------------------
@Pigolotti2013

The drastic assumptions of neutral theory favour mathematical tractability leading to analytical pre- dictions for ecological patterns such as species abundance distribution (see, e.g., Rosindell et al, 2011, and references therein). Such predictions, depending on few parameters, fit surprisingly well field data from tropical forests (see, e.g., Hubbell, 2001; Volkov et al., 2003, 2005). Due to its strong departure from traditional theoretical approaches, the neutral theory elicited a heated debate in community ecology about its validity and interpretation. A complication is that, in many cases, niche-based and neutral models yield similar fits of biodiversity patterns (see, e.g., Chave et al., 2002; McGill, 2003; Mouquet and Loreau, 2003; Tilman, 2004; McGill et al., 2006)

--------------------------------
@Stauffer1994

The divergence of characteristic times (in our case the fire lifetime) at a critical point it is called "critical slowing down". For example for a temperature only slight below the liquid gas critical temperature, the fluid is quite unsure whether it wants to be liquid or vapour, and thus take a lot of time to make its choice. 


--------------

@Sole2004

This occurs at the percolation value D ∗ ≈ 0.41 (i.e. when the proportion of avail- able sites is 1 − D ∗ ≈ 0.59) where it can be shown that the largest available patch cannot longer be connected and splits into many sub-patches. The same effect has been found in our system. In Fig. 5(a), we show the biodiversity decay linked to the presence of a percolation phenomenon in the spatially-explicit landscape. As we can see, there is a steady decay in the number of species present as D increases, with a rapid decay as D ∗ is approached. A different view of the decay can be obtained by plotting the Shannon entropy (Fig. 5(c)). We can see that it first increases as the best colonizers become less and less abundant (and the species-abundance distribution becomes more uniform). At D ∗ the effective loss of diversity reduces the value of the entropy which drops rapidly. As it was shown by the MF model, the average colonization rate also decreases (Fig. 5(b)). As D ∗ is approached, the average colonization rate goes down.


---------
@Stauffer1979 Scaling theory of percolation clusters

Thus percolation can be used as a guide to clustering phenomena at other phase transitions. Of course
there is no guarantee that other phase transitions will have the same cluster properties as the simple
percolation problem; but it seems plausible that a complete understanding of simple clusters (i.e.
percolation) is helpful for a better understanding of more complicated clusters (e.g. fluids, magnets).


----------
@Sornette2009

Several previous works [560, 774] have used a method to determine the
critical percolation concentration which uses a similar feedback procedure
(D. Stauffer, private communication): occupy a lattice with some probabil-
ity p and check if it percolates. If yes, decrease p by δp, if not increase p by
δp. Then, occupy the lattice again using the same sequence of random num-
bers and the new p. Then check and decrease δp by a factor 2, and change
p again by the new δp. Thus, in a square site problem (p c = 0.593...), one
may obtain the sequence (starting with p = 1/2, δp = 1/4): p = 1/2, no;
p = 3/4, yes; p = 5/8, yes; p = 9/16, no... This is an exponentially efficient
iteration to determine p c for one fixed configuration. The method has been
called “Hunting the Lion” by A. Aharony: you move in the direction where
you hear the lion roar, the lion being here the “infinite” percolating cluster.


------------
@Windus2007

Finding the critical exponents through steady state simulations is notoriously difficult
due to critical slowing down, finite-size effects, large fluctuations and the difficulties
that arise in finding the critical point. A much more effective method is that of time-
dependent simulations, which has proved to be a very efficient way of determining the
critical exponents and the critical point for models exhibiting absorbing phase transitions
[13]


------------
@Bregman2015 Species interactions regulate the collapse of biodiversity and ecosystem function in tropical forest fragments

Here, we use phylogenetic and functional trait data to test whether communities of two ecologically important guilds of tropical birds (frugivores and insectivores) are structured by species interactions in a fragmented Amazonian forest landscape. In both guilds, we found that forest patch size, quality, and degree of isolation influence the phylogenetic and functional trait structure of communities, with small, degraded, or isolated forest patches having an increased signature of competition (i.e., phylogenetic and functional trait overdispersion in relation to null models). These results suggest that local extinctions in the context of fragmentation are nonrandom, with a consistent bias toward more densely occupied regions of niche space. We conclude that the loss of biodiversity in fragmented landscapes is mediated by niche-based competitive interactions among species, with potentially far-reaching implications for key ecosystem processes, including seed dispersal and plant damage by phytophagous insects.


-------------
@Campos2013

The effects of habitat fragmentation and their implications for biodiversity is a central issue in conservation biology which still lacks an overall comprehension. There is not yet a clear consensus on how to quantify fragmentation even though it is quite common to couple the effects of habitat loss with habitat fragmentation on biodiversity. Here we address the spatial patterns of species distribution in fragmented landscapes, assuming a neutral community model. To build up the fragmented landscapes, we employ the fractional Brownian motion approach, which in turn permits us to tune the amount of habitat loss and degree of clumping of the landscape independently. The coupling between the neutral community model, here simulated by means of the coalescent method, and fractal neutral landscape models enables us to address how the species-area relationship changes as the spatial patterns of a landscape is varied. The species-area relationship is one of the most fundamental laws in ecology, considered as a central tool in conservation biology, and is used to predict species loss following habitat disturbances. Our simulation results indicate that the level of clumping has a major role in shaping the species-area relationship. For instance, more compact landscapes are more sensitive to the effects of habitat loss and speciation rate. Besides, the level of clumping determines the existence and extension of the power-law regime which is expected to hold at intermediate scales. The distributions of species abundance are strongly influenced by the degree of fragmentation. We also show that the first and second commonest species have approximately self-similar spatial distributions across scales, with the fractal dimensions of the support of the first and second commonest species being very robust to changes in the spatial patterns of the landscape.

-------------
@Fort2013

The classical competition niche theory (CCNT) based on the Lotka-Volterra competition equations and Hutchinson’s multidi- mensional niche [1,2]

Being far from equilibrium, fitting detailed snapshots of tropical forest data can only describe transitory configurations but cannot help to understand the mechanisms underlying the observed dynamics. The second problem of both CCNT and NTB is the reliance on mean-field approximations stemming from the ‘‘well-mixed’’ assumption that clearly fails when applied to sessile individuals like trees having largely local recruitment and potentially suffering from local interference competition for resources. Thus, neither CCNT nor NTB are suitable for the realistic, parsimonious modelling of tropical forest dynamics, nor can they predict the large array of metrics characterizing the spatio-temporal dynamics of these forest communities.


--------------

@Chisholm2009c   Linking dispersal, immigration and scale in the neutral theory of biodiversity

The immigration parameter m is the probability that a death in the local community is replaced by the offspring of an individual from outside the local community.

the mean dispersal distance d, based on seed-trap data for 81 tree species in the plot, is 39.5 m and the 95% confidence interval for d is (32.8, 46.7) (confidence interval estimated by bootstrap-resampling raw data in Muller-Landau 2001; see also Condit et al. 2002)

Although m is useful for analytical purposes, such as deriving the theoretical functional forms of SADs (e.g., Volkov et al. 2003), we propose that the mean dispersal distance d is a more suitable fundamental parameter because it is scale invariant, biologically meaningful and easily related to m analytically

The predictions of spatially implicit neutral theory appear robust to violations of the assumption that all species in a community have the same dispersal kernel: our main approximation (2) holds even when different species have different dispersal kernels and exert different propagule pressures, provided that the mean dispersal distance d is calculated as the average of the mean dispersal distances of different species (weighted by the speciesÕ relative abundances and propagule pressures; see Appendix S2). This is important, because the assumption that species have identical dispersal kernels is patently false

Our finding of robustness to interspecific variation in dispersal kernels provides some justification for a fundamental premise of neutral theory, which is that inter and intraspecific trait variation averages out statistically at large scales and is not important for predicting macroscopic patterns of diversity (Chave 2004, p250–251).

-----------------

@Etienne2011c The spatial limitations of current neutral models of biodiversity.

binned abundance distributions are largely uninformative unless the visual fit is bad, in which case the statistical fit will be bad too. Rank-abundance distributions contain the full abundance vector of a sample and are therefore more informative. This has been argued before (e.g. [35]), but nevertheless binned SADs are still widely used. Perhaps our results will help to make this point again.



-----------------
@Kefi2013 Early warning signals also precede non-catastrophic transitions

Due to statistical limitations, a trend in a single indicator is likely to be insufficient to determine the proximity to a transition in real ecosystems (Biggs et  al. 2009). Therefore combining several indicators to confidently detect an approach- ing transition (Guttal and Jayaprakash 2009, Drake and Griffen 2010, Dakos et al. 2011), or using both time series and spatial information if available (Rietkerk et  al. 2004, Kéfi et  al. 2007, 2011) can be a promising way forward. Also combining early warning signals with predictions of mechanistic models may be helpful (Lenton 2011).

-----------------
@Aschehoug2015 Diversity Increases Indirect Interactions, Attenuates the Intensity of Competition, and Promotes Coexistence

In a field experiment, we found that competition among native perennial plants in multispecies assemblages was far weaker than competition between those same species in pairwise arrangements and that indirect interactions appeared to weaken direct competitive effects

-----------------
@McGill2010a Towards a unification of unified theories of biodiversity

Presumably if we can find the minimally sufficient set of rules such that the in papyro or in silico analyses match the in situ (field-based real world) analyses to accurately reproduce the macroecological patterns of biodiversity, we will have achieved a useful description of rules governing nature. This is the central goal of this paper.


-----------------
@Grilli2017 Feasibility and coexistence of large ecological communities

A particularly interesting feature of random matrices is that the distribution of their eigenvalues (determining local stability) is universal26. This means that local stability depends on just a few, coarse-grained properties of the matrix (that is, the number of species and the first few moments of the distribution of interaction strengths) and not on the finer details (for example, the particular distribution of interaction strengths)


We obtain that the mean of interaction strengths sets the behaviour of feasibility with the number of species. If the mean is negative (for example, in case of competition or predation with limited efficiency), the larger the system is, the smaller is the set of conditions leading to coexistence, while for positive mean (for example, in the case of mutualism) the converse is true.

Here we have shown that the fraction of conditions compatible with coexistence is mainly determined by the number and the mean strength of interactions. In terms of network properties, the relevant quantity is the connectance, with other properties (for example, nestedness or degree distribution) having minimal effects. In particular, once the connectance and mean interaction strength are fixed, the matrices built using empirical mutualistic networks have feasibility domains very similar to that expected for the random case, as was also observed previously
in a similar contex.

The major role is played by corse-grained statistical properties of the interactions, such as connectance or the mean and variance of the interaction strengths.

Consistently with earlier results 7,8 , this fact establishes a relationship between niche overlap and the range of conditions that lead to coexistence: greater niche overlap means a more restricted parameter range allowing for coexistence, irrespective of the details of the interactions.