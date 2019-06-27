# Updating Rosetta Interaction Potentials with Topology Derived Amino-Acid Pairwise Statistics

## Motivation
* Rosetta structure prediction and modeling is very sophisticated. 
* Modeling membrane proteins, however, present unique challenges for modeling in Rosetta from modeling soluble protein domains. 
* An important aspect to Rosetta modeling strategy is used statistically-derived residue pair-wise interaction potentials. 
* These potentials are derived from statistical analysis of pairwise contacts in experimentally determined protein structures. 
* Interacting pairs are often determined from using distance cutoffs and constraints. 
* Determining contacts from distances can be misleading, or over parameterized due to use of knowledge-based distance cutoffs. However, more importantly, a distance based view of residue interactions often misses the topology of residue interaction networks in proteins. 
* With distance cutoffs can over-estimate residue interactions. Using a topology based approach can more accurately parse true structural contacts from protein structures. 
* Building interaction potentials from this point of view may significantly influence the ability for Rosetta to identify correct folds of membrane proteins.  
