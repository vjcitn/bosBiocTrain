bosBiocTrain
============

Training materials for Bioconductor courses in Longwood area, Boston MA

October 10-11 2013, Fenway Room, Inn at Longwood
432 Longwood Ave Boston MA 02115

Day 1
=====

 * Lecture 1: Pitfalls of genomic data analysis.
  Complexity, poor design, batch effects.  Vehicles for
  avoiding some of the pitfalls with Bioconductor.

 * Lab 1: Managing genomic annotation for humand and model organisms
   + OrganismDbi, select
   + GenomicRanges

 * Lab 2: Managing and using experimental archives
   + General principle: X[G, S] is selection of genomic features G
    and experimental samples S from archive X
   + Archive containers: ExpressionSet, SummarizedExperiment
   + Getting acquainted with some classic experiments
     - Expression arrays
     - Methylation arrays
     - Genotyping studies
     - NGS studies: RNA-seq, ChIP-seq
 * Lecture 2: Statistical concepts for genomic data analysis
   + Exploratory data analysis
     - distributions, density estimation
     - scatterplot matrices
     - PCA
     - distances, clustering, silhouette
     - Example: identifying batch effects
   + Hypothesis testing
     - Two-sample problem: parametric, nonparametric
     - regression/ANOVA
     - censored response
     - correlated response
   + Shrinkage concepts for high-dimensional data
   + Visualizations and reports: standard, "shiny", ReportingTools

 * Lab 3: Statistical explorations of genomic data: interfaces
for exploratory multivariate analysis, machine learning,
multiple comparisons, enumerating significantly distinctive
features, functional interpretation of feature sets

Day 2
=====
 * Lecture 3: R and Bioconductor for high-throughput computing

 * Lab 4: Case studies
   + Microarray differential expression
   + Differential methylation
   + RNA-seq, DESeq
   + ChIP-seq, bsseq, DiffBind
   + Integrative analyses
