# The Skeletal Cell Atlas

## About
Breakthroughs in improving scalability have led to the creation of organism-wide transcriptomic datasets, aiming to comprehensively profile the cell types and states within an organism throughout its lifecycle. To date however, the skeleton remains a majorly underrepresented organ system in organism-wide atlases. Considering that the skeleton is the central framework of the vertebrate body, the home of the hematopoietic niche and a central player in major metabolic and homeostatic processes, this presents a major deficit in current reference atlas projects. To address this issue, we integrated seven separate scRNA-seq datasets containing skeletal cells and their developmental precursors, generating an atlas of over 800,000 cells. This skeletal cell atlas describes cells across the mesenchymal lineage from the induction of the limb field to adult bone, encompassing 50 different cell states.

This repository contains all notebooks required to perform the analyses performed in the manuscript and generate all figures.

## Accessibility 
The atlas can be freely explored at http://www.skeletalcellatlas.org. The code used to create this webapp can be found [here](https://github.com/mbarzegary/skeletal-cell-atlas).
## Reproducibility
To avoid conflicts between package dependencies, it is strongly recommended to use our [Docker](https://hub.docker.com/r/gnasello/sc-env) environment.
