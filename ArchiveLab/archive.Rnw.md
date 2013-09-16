%\VignetteIndexEntry{test rnw vignette}
%\VignetteEngine{knitr::knitr}


<h2>Representing genome-scale data with R</h2>

For most practical purposes, ``genomes" are big data, and
even if we consider management of genomic sequence data to be
well under control, management of the refinement of reference
sequences and annotations at the species or individual level
rapidly involves us in substantial problems of complexity and
volume.  The internet is a unifying principle, and some approaches
to large-scale distributed computing (Amazon EC2, Google Compute)
confer some unification on the approaches that will often be taken
to deal with very large-scale genomic data processing.

We adopt the R programming language for the representation and
analysis of genome-scale data for many reasons.  Some are listed here.
<ul>
<li>Many statisticians implement new algorithms for inference in R
for the purpose of testing, refining, illustrating, and distributing new
methods.  We want access to these methods for our genomic research.
<li>R is freely distributed and redistributable, and is maintained
in a freely redistributable form on a number of important hardware and
software platforms.
<li>Contributions to the R data analysis environment foster interoperability
with other software tools including relational
databases, tools for network visualization and
inference, sequence management software,
and high-performance computation, so it is seldom the case
that use of R is an important restriction on data-analytic versatility.
</ul>

In any event, we'll consider some very basic issues of representation
here.


<h3>Introduction to S4 object-oriented programming</h3>

Object-oriented programming is a style of programming that emphasizes
the formalized
unification of heterogeneous data in "objects", defining relationships
among classes of objects, and promoting software designs in which
function behavior can vary according to the classes of objects on
which the function is evaluated.  In R's "S4" object-oriented methodology,
we use 
<ul>
<li><code>setClass( [classname], [representation], ...) </code> to
define a class of objects and the classes of entities from which instances
of the class are composed,
<li><code>setGeneric( [genericName], [interface], ...) </code> to
establish a generic function whose behavior will depend upon
the classes of objects on which it will be evaluated,
<li><code>setMethod( [methodName], [signature], ...) </code>  to
establish the method for a specific ordered combination of objects (the ``signature'')
on which the generic function is to be evaluated.
</ul>
For full details, see monographs by Chambers (Software for Data Analysis)
and Gentleman (R Programming for Bioinformatics).

<h3>Genomic sequence: BSgenome</h3>

In this section we use the \textit{Biostrings}
infrastructure to check the genomic location of
a yeast microarray probe.

First, we obtain the reference genomic sequence
for sacCer2 version of yeast.

```r
library(BSgenome.Scerevisiae.UCSC.sacCer2)
```

```
## Loading required package: BSgenome Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
## clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport,
## clusterMap, parApply, parCapply, parLapply, parLapplyLB, parRapply,
## parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
## xtabs
## 
## The following objects are masked from 'package:base':
## 
## anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
## duplicated, eval, Filter, Find, get, intersect, lapply, Map, mapply,
## match, mget, order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
## rbind, Reduce, rep.int, rownames, sapply, setdiff, sort, table, tapply,
## union, unique, unlist
## 
## Loading required package: IRanges Loading required package: GenomicRanges
## Loading required package: XVector Loading required package: Biostrings
```

```r
class(Scerevisiae)
```

```
## [1] "BSgenome"
## attr(,"package")
## [1] "BSgenome"
```

```r
Scerevisiae
```

```
## Yeast genome
## | 
## | organism: Saccharomyces cerevisiae (Yeast)
## | provider: UCSC
## | provider version: sacCer2
## | release date: June 2008
## | release name: SGD June 2008 sequence
## | 
## | sequences (see '?seqnames'):
## |   chrI     chrII    chrIII   chrIV    chrV     chrVI    chrVII   chrVIII
## |   chrIX    chrX     chrXI    chrXII   chrXIII  chrXIV   chrXV    chrXVI 
## |   chrM     2micron  
## | 
## | (use the '$' or '[[' operator to access a given sequence)
```


We focus attention on chrIV for now.

```r
c4 = Scerevisiae$chrIV
class(c4)
```

```
## [1] "DNAString"
## attr(,"package")
## [1] "Biostrings"
```

```r
c4
```

```
##   1531919-letter "DNAString" instance
## seq: ACACCACACCCACACCACACCCACACACACCACA...AAACATAAAATAAAGGTAGTAAGTAGCTTTTGG
```


Now we obtain the probe and annotation data for
the Affy yeast2 array.

```r
library(yeast2.db)
```

```
## Loading required package: AnnotationDbi Loading required package: Biobase
## Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
## 
## Attaching package: 'AnnotationDbi'
## 
## The following object is masked from 'package:BSgenome':
## 
## species
## 
## Loading required package: org.Sc.sgd.db Loading required package: DBI
```

```r
library(yeast2probe)
yeast2probe[1:3, ]
```

```
##                    sequence   x   y Probe.Set.Name
## 1 GAAAGTTTCAGTGCACGTCTTCAAA 380 257     1769438_at
## 2 GTATATTTCTAATCTTCCTCTTCAT  28 327     1769438_at
## 3 ATATCAAACCGCGTACTTCGTGACT 188  19     1769438_at
##   Probe.Interrogation.Position Target.Strandedness
## 1                         1117           Antisense
## 2                         1170           Antisense
## 3                         1240           Antisense
```


We'll pick one probe set, and then the sequence of one probe.

```r
ypick = yeast2probe[yeast2probe$Probe.Set.Name == "1769311_at", ]
dim(ypick)
```

```
## [1] 11  6
```

```r
ypick[1:3, ]
```

```
##                       sequence   x   y Probe.Set.Name
## 1959 ATGAGCACTATGTTTTCTGTTGGAT 486  39     1769311_at
## 1960 GTTTTCTGTTGGATTTGGCTCATAC 154 321     1769311_at
## 1961 TTGGCTCATACTTGGCATCTGGGAA  20 493     1769311_at
##      Probe.Interrogation.Position Target.Strandedness
## 1959                          100           Antisense
## 1960                          111           Antisense
## 1961                          125           Antisense
```

```r
a = "ATGAGCACTATGTTTTCTGTTGGAT"
ra = reverseComplement(DNAString(a))
```


Now we use the simple lookup of Biostrings.

```r
matchPattern(ra, c4)
```

```
##   Views on a 1531919-letter DNAString subject
## subject: ACACCACACCCACACCACACCCACACACACCA...ACATAAAATAAAGGTAGTAAGTAGCTTTTGG
## views:
##      start    end width
## [1] 174478 174502    25 [ATCCAACAGAAAACATAGTGCTCAT]
```

```r
get("1769311_at", yeast2CHRLOC)
```

```
##       4 
## -174232
```

```r
get("1769311_at", yeast2CHRLOCEND)
```

```
##       4 
## -174588
```


Exercise.  Discuss how to identify probes harboring SNP in the
affy u133plus2 array.


<h3>ExpressionSet and self-description</h3>

The <code>ExpressionSet</code> class unifies information on a microarray experiment.

```r
library(ALL)
data(ALL)
getClass(class(ALL))
```

```
## Class "ExpressionSet" [package "Biobase"]
## 
## Slots:
##                                                                
## Name:      experimentData          assayData          phenoData
## Class:              MIAME          AssayData AnnotatedDataFrame
##                                                                
## Name:         featureData         annotation       protocolData
## Class: AnnotatedDataFrame          character AnnotatedDataFrame
##                          
## Name:   .__classVersion__
## Class:           Versions
## 
## Extends: 
## Class "eSet", directly
## Class "VersionedBiobase", by class "eSet", distance 2
## Class "Versioned", by class "eSet", distance 3
```

```r
experimentData(ALL)
```

```
## Experiment data
##   Experimenter name: Chiaretti et al. 
##   Laboratory: Department of Medical Oncology, Dana-Farber Cancer Institute, Department of Medicine, Brigham and Women's Hospital, Harvard Medical School, Boston, MA 02115, USA. 
##   Contact information:  
##   Title: Gene expression profile of adult T-cell acute lymphocytic leukemia identifies distinct subsets of patients with different response to therapy and survival. 
##   URL:  
##   PMIDs: 14684422 16243790 
## 
##   Abstract: A 187 word abstract is available. Use 'abstract' method.
```


A very high-level aspect of self-description is the facility for
binding information on the publication of expression data to
the data object itself.  The <code>abstract</code> method gives
even more detail.

Fundamental operations on <code>ExpressionSet</code> instances are
<ul>
<li><code>[</code> : subsetting, so that <code>X[G, S]</code> restricts
the object to genes or probes enumerated in <code>G</code> and samples
enumerated in <code>S</code>,
<li><code>exprs()</code>: return G &times; N matrix of expression values
<li><code>pData()</code>: return N &times; R data.frame instance with
sample-level data
</ul>

To enumerate additional formal methods, try
<code>
showMethods(classes="eSet", inherited=FALSE, showEmpty=FALSE, 
   where="package:Biobase")
</code>

<h3>SummarizedExperiment, VCF</h3>

Microarrays were used primarily with named probes.  You would find
differential intensity for <code>1007\_s\_at</code> and look that token
up.  With short read sequencing or very dense arrays, it 
is more relevant to have access
to the genomic location of the read or probe for interpretation.  The
<code>SummarizedExperiment</code> class addresses this concern.

A large-scale illustration of this class is in the \textit{dsQTL}
package, but it is very large so I do not require that you download it.

```r
library(dsQTL)
data(DSQ_17)
DSQ_17
rowData(DSQ_17)[1:3]
```


You can use \verb+example("SummarizedExperiment-class")+ to
get a working minimal example.

We use range-based operations to subset features; column
subscripting still selects samples.


The <code>VCF</code> class extends <code>SummarizedExperiment</code>
to manage information on variant calls on a number of samples.

```r
library(VariantAnnotation)
```

```
## Loading required package: Rsamtools
## 
## Attaching package: 'VariantAnnotation'
## 
## The following object is masked from 'package:base':
## 
## tabulate
```

```r
fl <- system.file("extdata", "structural.vcf", package = "VariantAnnotation")
vcf <- readVcf(fl, genome = "hg19")
vcf
```

```
## class: CollapsedVCF 
## dim: 7 1 
## rowData(vcf):
##   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
## info(vcf):
##   DataFrame with 10 columns: BKPTID, CIEND, CIPOS, END, HOMLEN, HOMSEQ,...
## info(header(vcf)):
##              Number Type    Description                                  
##    BKPTID    .      String  ID of the assembled alternate allele in th...
##    CIEND     2      Integer Confidence interval around END for impreci...
##    CIPOS     2      Integer Confidence interval around POS for impreci...
##    END       1      Integer End position of the variant described in t...
##    HOMLEN    .      Integer Length of base pair identical micro-homolo...
##    HOMSEQ    .      String  Sequence of base pair identical micro-homo...
##    IMPRECISE 0      Flag    Imprecise structural variation               
##    MEINFO    4      String  Mobile element info of the form NAME,START...
##    SVLEN     .      Integer Difference in length between REF and ALT a...
##    SVTYPE    1      String  Type of structural variant                   
## geno(vcf):
##   SimpleList of length 4: GT, GQ, CN, CNQ
## geno(header(vcf)):
##        Number Type    Description                                  
##    GT  1      String  Genotype                                     
##    GQ  1      Float   Genotype quality                             
##    CN  1      Integer Copy number genotype for imprecise events    
##    CNQ 1      Float   Copy number genotype quality for imprecise...
```


<h3>Management of information on HTS experiments</h3>


We'll get into details on this later.  For now, consider what aspects
of the structures you've encountered so far appear relevant to the
structure and management of HTS data.

