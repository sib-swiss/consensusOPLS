## ConsensusOPLS: An R package for Multi-Block Data Fusion

This is the github repository for the *ConsensusOPLS* package - an R
implementation and improvement of the Consensus Orthogonal Partial Least
Squares method, which was previously developed in Matlab
(https://gitlab.unige.ch/Julien.Boccard/consensusopls). The latest release of
the package is available on
[CRAN](https://cran.r-project.org/package=ConsensusOPLS).

Merging data from multiple sources is a relevant approach for comprehensively
evaluating complex systems. However, the inherent problems encountered when
analyzing single tables are amplified with the generation of multi-block
datasets, and finding the relationships between data layers of increasing
complexity constitutes a challenging task. For that purpose, a generic
methodology is proposed by combining the strength of established data analysis
strategies, i.e. multi-block approaches and the OPLS framework to provide an
efficient tool for the fusion of data obtained from multiple sources. The
package enables quick and efficient construction of the consensus OPLS
model for any horizontal multi-block data structures (observation-based
matching). Moreover, it offers an interesting range of metrics and graphics to
help to determine the optimal number of components and check the validity of
the model through permutation tests. Interpretation tools include score and
loading plots, Variable Importance in Projection (VIP), functionality predict
for SHAP computing and performance coefficients such as R2, Q2 and DQ2
coefficients.

This package is based on the *K-OPLS* algorithm, developed by Max Bylesjo,
University of Umea, Judy Fonville and Mattias Rantalainen, Imperial College.
Copyright (c) 2007-2010 Max Bylesjo, Judy Fonville and Mattias Rantalainen.
This *K-OPLS* code has been extended and adapted under the terms of the GNU
General Public License version 2 as published by the Free Software Foundation.

## Installation

```
# Installation via CRAN
install.packages("ConsensusOPLS")

# Installation via GitHub
devtools::install_github("sib-swiss/ConsensusOPLS")
```

## Bug report

Users can report bugs and request new features through the GitHub issue system
or by directly emailing the package maintainer for continuous improvement and
support.
	
## License

This package is licensed under GPL (>=3).

## Authors

Contributors: Céline Bougel, Julien Boccard, Florence Mehl, Marie Tremblay-Franco,
Mark Ibberson, Van Du T. Tran

## Acknowledgments

This project was developed in collaboration between:

- Institut national de recherche pour l'agriculture, l'alimentation et
l'environnement (INRAE), within the Food Toxicology unit (ToxAlim - Axiom
technical platform) - funded by ANR MetaboHUB (MTH), the National Metabolomics
and Fluxomics Infrastructure version 2.0.

- SIB Swiss Institute of Bioinformatics - Vital-IT.

