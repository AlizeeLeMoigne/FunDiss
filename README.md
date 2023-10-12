# FunDiss
The code descibed in this project is linked with the manuscript: A global eco-evolutionary pattern of bacterial functional diversity
Manuscript Authors: Alizée Le Moigne, Adrian-Stefan Andrei and Jakob Pernthaler
Scripts Author: Alizée Le Moigne

The datasets to be used with the script are publicly available at: https://doi.org/10.6084/m9.figshare.c.6734349

Description of the project: Microbial communities are at the heart of most biogeochemical processes on Earth, empowered by the astounding diversity of their functional traits. However, a comprehensive framework to address the functional responses of microbial communities
to environmental change is still lacking.  Through the integration of theoretical ecology with large-scale genome-resolved metagenomics, we have generated a systematic overview of the diversity of bacterial functional profiles. Upon the evaluation of major categories of
functional traits, we found higher dissimilarity among genomes in traits related to the interaction with the environment than in traits tied to basic functions like DNA replication and translation. Additionally, this distinct pattern across functions was modulated by the
microbiome type, highlighting its potential to identify fundamental mechanisms governing functional diversity across habitats, ecosystems, and communities. Interestingly, aquatic microbiomes exhibited higher functional diversity than the other microbiomes. The functional
diversity pattern described in this study offers a universal framework for understanding and model ecosystem-wide dynamics of microbial communities.

## 1. 1Fundiss_microcosms.R
This R script was used to calculate the functional dissimilarity at various biodiversity scales (alpha, beta and gamma) of bacterial communities grown at identical environmental conditions. The bacteria were collected from 5m depth of Lake Zurich (Switzerland) and grown in
the laboratory. After two weeks, the DNA of these bacterial communities was collected and sequenced (metagenomics). The plots corresponds to Fig 1 Panel A of the manuscript. The statistical analysis of the functional dissimilarity of the communities from the microcosms 
are available in this script.

## 2. 2GTDB_null.R
This script was used to infer the functional dissimilarity of publicly available bacterial genomes from the Genome Taxonomy DataBase (GTDB) originating from various environments. The plot correspond to Fig 1 Panel B of the manuscript. The statistical analysis of the
publicly available genomes (GTDB) is available in this script.

## 3. 3Env_specificity.R
This script was used to compare the functional dissimilarity of bacterial genomes originating from different environments/habitat. The genomes were extracted from the GTDB. The plots corresponds to Fig 2 of the manuscript. The statistical analysis , i.e. comparison of
functional dissimilarity between different environments is available in this script.

## 4. modified_Rao.R
This function was modified from the rao.diversity() available in the R package SYNCSA. The modification consists in replacing the Gower distance (Euclidean) with the Bray-Curtis dissimilarity in the calculation of the functional dissimilarity (dissimilarity of the species
traits matrices), Bray-Curtis being more adapted to sparse matrices.
