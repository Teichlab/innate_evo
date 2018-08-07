Gene expression variability across cells and species shapes innate immunity
===========================================================================

This directory contains data files and scripts related to our work: 
"Gene expression variability across cells and species shapes innate immunity".

In this work, we have studied the relationship between transcriptional divergence
between species in response to immune challenges and expression variability between
individual cells.  We used two cross-species systems to study this: dermal fibroblasts
and bone marrow-derived mononuclear phagocytes.

We further study the relationship between these variability and divergence and other
characteristics - including: 
1. promoter architecture (presence or absence of elements such as CpG islands or TATA-boxes),
2. gene functions, and 
3. other evolutionary and cellular characteristics (such as coding sequence evolution, number of cellular interactions).

The characteristics used are always related to one reference species: human - in the case of
fibroblasts, and mouse - in the case of phagocytes. For example, when studying the relationship
between cross-species divergence and promoter architecture, we use human and mouse promoters,
for fibroblasts and phagocytes, respectively.

The scripts in this directory allow reproducing the statistical analysis and main plots.
We also provide a script detailing how we compute cell-to-cell variability using the DM method.
Larger files and raw data (such as fastq files and count matrices) can be found in our
FTP website and the relevant accessions in ArrayExpress (see manuscript for details).

Most groups were compared using empirical p-values. Similar results are obtained using either
the Mann-Whitney U test or an alternative bootstrap method. The code can be changed to output
p-values by any method by changing a flag in common_functions.R
