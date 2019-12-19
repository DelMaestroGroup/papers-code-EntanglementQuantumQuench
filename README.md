<!-- [![Paper](https://img.shields.io/badge/paper-arXiv%3A1905.03312-B31B1B.svg)](https://arxiv.org/abs/1905.03312) -->

# Universal Entanglement Growth After a Quantum Quench
Adrian Del Maestro, Hatem Barghathi and Bernd Rosenow

### Abstract
The time evolution of an initial quantum state after a sudden change of interaction strength  leads to an asymptotic steady state, whose local properties are governed by the buildup of entanglement between spatial subregions of the system. Little is known about the simultaneous evolution of entanglement between groups of particles, which is based on n-point correlation functions.  We  analyze fermions after an interaction quantum quench in one spatial dimension and demonstrate that the steady state entropy density accumulated via entanglement dynamics is equivalent under either a spatial or particle bipartition.  Building on this connection, we  find that experimentally accessible density-density correlations can be employed to construct a  diagonal ensemble density matrix of  non-interacting bosons to compute the von Neumann entropy density, complimentary to measurements of  Renyi entropies. Our results highlight the universality of the dynamical transmutation of entanglement to thermodynamic entropy under time evolution that underlies our current framework of quantum statistical mechanics.

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project was generated via exact diagonalization.  Everything included in the [data](https://github.com/DelMaestroGroup/papers-code-EntanglementQuantumQuench/tree/master/data) directory was generated via:

* [ED for time evolution and entanglement entropy](https://github.com/DelMaestroGroup/tVDiagonalizeTimeEvaluationQuench/tree/TranslationalSymmetricInitialState_IntFermionBasis)
* [LMFIT for finite size scaling analysis](https://lmfit.github.io/lmfit-py/)

### Support
The creation of these materials was supported in part by the National Science Foundation under Award No. DMR-1553991.

[<img width="100px" src="https://www.nsf.gov/images/logos/NSF_4-Color_bitmap_Logo.png">](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1553991)

### Figures

#### Figure 01: Entanglement evolution after a quantum quench
<img src="https://github.com/DelMaestroGroup/papers-code-EntanglementQuantumQuench/blob/master/figures/method_flowchart_nobox.svg" width="600px">

#### Figure 02: Exact diagonalization results for entanglemnent
<img src="https://github.com/DelMaestroGroup/papers-code-EntanglementQuantumQuench/blob/master/figures/DeltaS_vs_t_ED.svg" width="600px">

#### Figure 03: Universality of entanglement
<img src="https://github.com/DelMaestroGroup/papers-code-EntanglementQuantumQuench/blob/master/figures/asymptotic_entropy_vs_invN_LL_prediction.svg" width="600px">

#### Figure 04: Entanglement entropy from the bosonic diagonal ensemble
<img src="https://github.com/DelMaestroGroup/papers-code-EntanglementQuantumQuench/blob/master/figures/nobdm_SqnLL_N12_EE.svg" width="600px">

