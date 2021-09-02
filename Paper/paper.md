---
title: 'A coupled model for the accumulation of cosmogenic radionuclides during shore platform development'
tags:
  - Geomorphology
  - Rock Coasts
  - Cosmogenic Radionculdies
  - C++

authors:
  - name: Martin D. Hurst^[co-first author]^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-9822-076X
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Hironori Matsumoto^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-4637-7781
    affiliation: 2
  - name: Jennifer R. Shadrick^[co-first author]
    affiliation: 3
  - name: Mark E. Dickson^[co-first author]
    orcid: 0000-0002-6031-4170
    affiliation: 4

affiliations:
 - name: School of Geographical and Earth Sciences, University of Glasgow, Glasgow, Scotland, UK.
   index: 1
 - name: Scripps Institution of Oceanography, University of California San Diego, La Jolla, CA, US.
   index: 2
 - name: Earth Science and Engineering, Imperial College London, London, UK.
   index: 3
 - name: School of Environment, Univesity of Auckland, Auckland, NZ.
   index: 4
date: 28th July 2021
bibliography: paper.bib

---

# Summary

Rock coasts evolve relatively slowly and oftentimes episodically, and because they are erosive in nature there is little evidence of their former state. 
These two factors compound to make it difficult to constrain the nature of processes that dictate their evolution, and associated rates of change across appropriately long timescales (centuries to millennia).
Measuring how rapidly rock coasts have evolved in the past is important to understand how they will react to future environmental change.
Recent developments in numerical modelling of rock coast weathering and erosion processes have shed new light on the fundamental controls on rock coast evolution at these timescales. 
In parallel, measurement of the concentration of rare cosmogenic radionuclides (CRNs) that accumulate in rocks at the coast is providing a new empirical basis for understanding the development of rock coasts over similarly long timescales. 
The accumulation of CRNs occurs most rapidly when rocks are near to the Earth surface, and thus the rates and processes by which rocks are unveiled at the coast is a first order control on the amount of CRNs that will be found in rock samples at the coast. 
Thus, numerical models that account for the morphological development of the shore platform coupled with modelled CRN production are vital if measured CRN concentrations are to reveal the style and pace of rock coast morphodynamics. 
Towards achieving this end, a coupled model of rock coast morphodynamics and CRN accumulation is presented here.

# Statement of need

Rock coasts make up a substantial proportion of the global coastline. 
Rock decay and wave-driven erosion combine to drive the landward retreat of bedrock, typically resulting in cliffed coasts fronted by shore platforms [@Kennedy2014].
The products of erosion (regolith/debris/beach material) are gradually destroyed or removed in this erosion-limited environment. 
Rates of erosion by cliff retreat and shore platform lowering are typically considered to be slow.
A recent compilation of global cliff retreat rates identified that rates vary from a few mm yr\textsuperscript{-1} to 10s cm yr\textsuperscript{-1}, with the subtstrate lithology exerting a dominant control on retreat rates [@Premaillon2018]. 
Shore platform erosion rates are typically on the order of mm yr\textsuperscript{-1} [@Stephenson2012; @Cullen2018; @Swirad2019]. Measurement of cliff retreat and shore platform erosion at rock coasts is temporally limited to at most a few decades. 
Therefore, reconstruction of the long-term (centennial-millennial) development, prior to observational records has until recently, been inferential at best, based on conceptual understanding of processes extrapolated using exploratory numerical models [@Trenhaile2000a; @Ashton2011; @Matsumoto2016a].

Rocks exposed at the Earth surface harbour an archive of cosmogenic radionuclides (CRNs) that reflect how long they have been uncovered or how rapidly they were unveiled to the surface [@Dunai2010; @Gosse2001]. 
The constant bombardment of the Earth by cosmic radiation causes nuclear reactions that produce CRNs when high energy cosmic ray particles collide with matter down to a few metres below the Earth surface. 
Measurement of CRN concentrations in bedrock and sediment samples are routinely used in thefield of Geosciences to constrain the age of bedrock outcrops or sedimentary deposits, or to estimate rates of erosion at a local or catchment scale [e.g. @Granger2013a]. 
Thus, concentrations of CRNs in shore platforms reflects the style and rate of morphological development of the coast [@Regard2012; @Hurst2017]. 

Several recent studies have measured CRN concentrations in rock samples collected from shore platforms in order to estimate long-term rates of cliff retreat or to demonstrate the antiquity of shore platforms at the coast [@Regard2012; @Choi2012; @Rogers2012; @Hurst2016; @Raimbault2018; @Swirad2020].
However, the distribution of CRN concentrations across shore platforms is dependent on the pace and style of morphological changes, and thus we require morphodynamic models that represent the myriad of processes influencing shore platform development coupled to predictions of CRN concentration distributions [@Hurst2017].
These previous efforts had represented the morphological development of rock coasts in a highly simplistic way, assuming the shoreward translation of a constant shore platform geometry, with the elevation of the cliff-platform junction tracking changes in mean relative sea level. 
However, there is little basis for this notion of equilibrium in rock coast development [@Dickson2013]. 
Numerical models of shore platform development suggests that shore platforms tend to widen through time [@Trenhaile2000a; Matsumoto2016a], with cliff retreat rates declining through time as wave energy is dissipated across the widening shore platform. 
The exception is when relative sea level rise results continues to facilitate wave energy reaching landward parts of the shore platform, such that cliff retreat continues in some way proportional to the rate relative sea level rise [@Ashton2011], and an equilibrium condition can emerge.

Morphodynamic models of rock coast development can exhibit a range of behaviours depending on intrinsic and extrinsic conditions including tidal range, relative sea level change, lithological properties and wave energy delivery [@Trenhaile2014; @Matsumoto2016a; @Matsumoto2018]. However, it is not well known how these behaviours impact upon the expected distribution of CRNs stored in the coastal bedrock.
Trenhaile suggested that “while [CRN] has the potential to revolutionize our understanding of the evolution of rock coasts, the accuracy of the results is dependent on the validity of the conceptual and mathematical model assumptions that have to be made” [-@Trenhaile2018, p. 80]. 
He goes on to highlight several key research challenges relating to processes that are poorly represented in morphodynamic models (such as weathering processes), or absent entirely (such as abrasion, block detatchment/quarrying and biological processes). 
This highlights a critical circular problem, that CRNs can reveal rates of rock coast development but requires a model of rock coast evolution that faithfully represents the dominant processes driving change, while measurements of CRNs has been suggested could improve our process understanding. 
Nevertheless, in order to progress efforts to quantify long-term rock coast development, a coupled modelling framework for the morphodynamic development of rock coasts and accumulation of cosmogenic radionuclides is required. 

# RPM-CRN

This paper contributes a coupled model combining rock coast and shore platform development [@Matsumoto2016a] and CRN production at rock coasts [@Hurst2017] in order to provide a framework for interpreting the history of coastal evolution from measured CRN concentrations. RPM-CRN is being used to intepret the topography and CRN concentrations at various rock coast settings globally through multi-objective optimisation to these field data [@Shadrick2021]. 

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Acknowledgements

This modelling effort has benefitted from discussions with Michael A. Ellis, Wayne J Stephenson and Larissa A Naylor. We acknowledge financial support from the Marine Alliance for Science and Technology Scotland (MASTS).

# References