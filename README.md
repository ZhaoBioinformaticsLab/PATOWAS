# PATOWAS
PATOWAS: A Pipeline for Analyzing Trait through 'Ome'- wide Association Studies


README

Updated on 06/01/2018


BACKGROUND

We proposed a novel linear mixed model (LMM) that is capable of accounting for not only the additive effects of biological markers but also the interaction effects of marker-by-marker pairs in multiple ‘ome’ (genome, transcriptome, proteome and metabolome)-wide association studies. 

Based on this LMM, we developed a pipeline, namely PATOWAS - A Pipeline for Analyzing Trait through 'Ome'- wide Association Studies, for analyzing traits through ome-wide association studies.  PATOWAS can therefore be used to study transcriptome- and metabolome-wide associations in addition to genome-wide associations. Such a technology is also called integrative ‘omics’.

The PATOWAS web server is publically available at http://bioinfo.noble.org/PATOWAS/.

The full detail about PATOWAS was described in Wenchao Zhang, Xinbin Dai, Shizhong Xu* (shizhong.xu@ucr.edu) and Patrick X. Zhao* (pzhao@noble.org). 2D association mapping and integrative omics analysis provides systems biology view in trait analysis, under revision.


PATOWAS SOURCE CODE

This package contains the source code of the following function modules in PATOWAS.

1. Module Description
PATOWAS includes four function calculation modules: a) kinship matrix calculation (km_calc), b) variance component analysis (vc_anal), c) ‘omics’-wide p-value scanning (1D scanning) for main additive effect (ps_main), and d) ‘omics’-wide p-value scanning (2D scanning) for interaction effect (ps_inter).

2. Function Description
PATOWAS was developed to analyze phenotypic trait through omics wide association studies, which need one of data matrix file such as 1) additive genotypic data, 2) gene expression data, or 3) metabolite abundance data as the input to calculate the 2 kinship matrices: ka, kaa, and furthermore need the trait related phenotypic data file to estimate the variance component ratios and p-value scanning for main additive and interaction effects.

3. Development Description
All the function modules were developed in C++ in a cross-open source IDE Code:: Blocks and have been successfully compiled and tested in Linux system. Users may download, modify and compile the source code, and further build their customized pipelines.
