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
```{r lkg}
library(BSgenome.Scerevisiae.UCSC.sacCer2)
class(Scerevisiae)
Scerevisiae
```

We focus attention on chrIV for now.
```{r doc4}
c4 = Scerevisiae$chrIV
class(c4)
c4
```

Now we obtain the probe and annotation data for
the Affy yeast2 array.
```{r getp}
library(yeast2.db)
library(yeast2probe)
yeast2probe[1:3,]
```

We'll pick one probe set, and then the sequence of one probe.
```{r getps}
ypick = yeast2probe[yeast2probe$Probe.Set.Name=="1769311_at",]
dim(ypick)
ypick[1:3,]
a = "ATGAGCACTATGTTTTCTGTTGGAT"
ra = reverseComplement(DNAString(a))
```

Now we use the simple lookup of Biostrings.
```{r getm}
matchPattern(ra, c4)
get("1769311_at", yeast2CHRLOC)
get("1769311_at", yeast2CHRLOCEND)
```

Exercise.  Discuss how to identify probes harboring SNP in the
affy u133plus2 array.


<h3>ExpressionSet and self-description</h3>

The <code>ExpressionSet</code> class unifies information on a microarray experiment.
```{r lkall}
library(ALL)
data(ALL)
getClass(class(ALL))
experimentData(ALL)
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
```{r lkdsq,eval=FALSE}
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
```{r lkv}
library(VariantAnnotation)
fl <- system.file("extdata", "structural.vcf", package="VariantAnnotation")
vcf <- readVcf(fl, genome="hg19")
vcf
```

<h3>Management of information on HTS experiments</h3>


We'll get into details on this later.  For now, consider what aspects
of the structures you've encountered so far appear relevant to the
structure and management of HTS data.

