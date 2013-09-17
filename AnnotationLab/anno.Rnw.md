%\VignetteIndexEntry{test rnw vignette}
%\VignetteEngine{knitr::knitr}

<h2>Annotation: genes, etc.</h2>
<h3>Prologue: programming with genomic coordinates</h3>

We would like to be able to consistently describe the
positions on the genome of various genomic and genetic
features.  The assembly of the human genome reference known
as GRCh37 (build 37 of the Genome Reference Consortium) is
also known as hg19 in nomenclature established in the UCSC
Genome Bioinformatics Group.  When referring to a genomic position,
it is essential to know which reference has been used to
define the address.  This is often not done: metadata defining
the context of positions is not concretely available, and
implicit understanding is assumed, necessitating checking
in case uncertainty is present.  It is good to be able to check:
Bioconductor SNPlocs packages have a long history and can
be used to verify, for SNP location assertions, which
genome build is in play.  However, it is better to use
objects that are self-describing and explicitly identify
the address basis.

Perhaps the most basic resource for working with
genomic content is the reference sequence itself.  We will
postpone this and work with higher-level concepts such
as organisms, genes, transcripts.

<h3>Genes and functional annotation, OrganismDbi</h3>

<h4>Overview</h4>

The OrganismDbi infrastructure defines relationships between
annotation resources useful in analyzing annotations for a specific
organism.


```r
library(Homo.sapiens)
```

```
## Loading required package: AnnotationDbi Loading required package:
## BiocGenerics Loading required package: parallel
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
## Loading required package: Biobase Welcome to Bioconductor
## 
## Vignettes contain introductory material; view with 'browseVignettes()'. To
## cite Bioconductor, see 'citation("Biobase")', and for packages
## 'citation("pkgname")'.
## 
## Loading required package: OrganismDbi Loading required package:
## GenomicFeatures Loading required package: IRanges Loading required
## package: GenomicRanges Loading required package: XVector Loading required
## package: GO.db Loading required package: DBI
## 
## Loading required package: org.Hs.eg.db
## 
## Loading required package: TxDb.Hsapiens.UCSC.hg19.knownGene
```

```r
Homo.sapiens
```

```
## class: OrganismDb 
## Annotation resources:
## [1] "GO.db"                             "org.Hs.eg.db"                     
## [3] "TxDb.Hsapiens.UCSC.hg19.knownGene"
## Annotation relationships:
##      xDbs           yDbs                                xKeys     
## [1,] "GO.db"        "org.Hs.eg.db"                      "GOID"    
## [2,] "org.Hs.eg.db" "TxDb.Hsapiens.UCSC.hg19.knownGene" "ENTREZID"
##      yKeys   
## [1,] "GO"    
## [2,] "GENEID"
## For more details, please see the show methods for the component objects listed above.
```


We see that the critical components of Homo.sapiens annotation
are 
org.Hs.eg.db.  
TxDb.Hsapiens...,  and
GO.db.
These address, respectively,
institutional cataloging and nomenclature
for genes, structural addressing of transcripts
and gene models, and functional, process, and anatomic classification
of genes and gene products.

Considerable metadata is provided for each resource.  For example,

```r
GO.db
```

```
## GODb object:
## | GOSOURCENAME: Gene Ontology
## | GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
## | GOSOURCEDATE: 20130302
## | Db type: GODb
## | package: AnnotationDbi
## | DBSCHEMA: GO_DB
## | GOEGSOURCEDATE: 2013-Mar5
## | GOEGSOURCENAME: Entrez Gene
## | GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
## | DBSCHEMAVERSION: 2.1
```


<h4>Lookups</h4>

A general task is taking a set of keys and finding their correspondences
in other resources.
The `columns' define names for
attributes of genomic elements that can be
used as types or tokens of keys or types of values to be requested.


```r
sort(columns(Homo.sapiens))
```

```
##  [1] "ACCNUM"       "ALIAS"        "CDSCHROM"     "CDSEND"      
##  [5] "CDSID"        "CDSNAME"      "CDSSTART"     "CDSSTRAND"   
##  [9] "CHR"          "CHRLOC"       "CHRLOCEND"    "DEFINITION"  
## [13] "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
## [17] "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "EXONCHROM"   
## [21] "EXONEND"      "EXONID"       "EXONNAME"     "EXONRANK"    
## [25] "EXONSTART"    "EXONSTRAND"   "GENEID"       "GENENAME"    
## [29] "GO"           "GOALL"        "GOID"         "IPI"         
## [33] "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL" 
## [37] "PATH"         "PFAM"         "PMID"         "PROSITE"     
## [41] "REFSEQ"       "SYMBOL"       "TERM"         "TXCHROM"     
## [45] "TXEND"        "TXID"         "TXNAME"       "TXSTART"     
## [49] "TXSTRAND"     "UCSCKG"       "UNIGENE"      "UNIPROT"
```

```r
s1 = AnnotationDbi::select(Homo.sapiens, keytype = "SYMBOL", keys = c("ABL1", 
    "BRCA2"), columns = c("UNIGENE", "ENTREZID", "UNIPROT", "ENSEMBL", "ENSEMBLPROT", 
    "ENSEMBLTRANS"))
```

```
## Warning: You have selected the following columns that can have a many to
## one relationship with the primary key: UNIGENE, ENSEMBL, ENSEMBLPROT,
## ENSEMBLTRANS, UNIPROT . Because you have selected more than a few such
## columns there is a risk that this selection may balloon up into a very
## large result as the number of rows returned multiplies accordingly. To
## experience smaller/more manageable results and faster retrieval times, you
## might want to consider selecting these columns separately. Warning:
## 'select' resulted in 1:many mapping between keys and return rows Warning:
## 'select' resulted in 1:many mapping between keys and return rows
```

```r
dim(s1)
```

```
## [1] 48  7
```

```r
s1[!duplicated(s1$SYMBOL), ]
```

```
##    SYMBOL ENTREZID   UNIGENE         ENSEMBL     ENSEMBLPROT
## 1    ABL1       25 Hs.431048 ENSG00000097007 ENSP00000361423
## 19  BRCA2      675  Hs.34012 ENSG00000139618 ENSP00000369497
##       ENSEMBLTRANS UNIPROT
## 1  ENST00000372348  P00519
## 19 ENST00000380152  P51587
```


Entrez gene identifiers are simple integers.  For working with
homology annotation resources it is useful to be acquainted with
Ensembl gene, transcript, and protein identifiers.

<h3>Gene models and genomic coordinates: TranscriptDb</h3>

To acquire coordinates of transcripts and exons, we use
a TranscriptDb instance.
The query is posed in terms of the target resource (here the hg19
UCSC transcript set), a <code>vals</code> specification (we give a binding
to <code>gene\_id}) and what we want back, in \texttt{col</code>.  We will get
a <code>GRanges</code> instance with metadata.

```r
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
tsel = transcripts(txdb, vals = list(gene_id = unique(s1$ENTREZID)), col = c("gene_id", 
    "tx_id", "tx_name"))
tsel
```

```
## GRanges with 4 ranges and 3 metadata columns:
##       seqnames                 ranges strand |         gene_id     tx_id
##          <Rle>              <IRanges>  <Rle> | <CharacterList> <integer>
##   [1]    chr13 [ 32889617,  32907524]      + |             675     48293
##   [2]    chr13 [ 32889617,  32973809]      + |             675     48294
##   [3]     chr9 [133589268, 133763062]      + |              25     34764
##   [4]     chr9 [133710831, 133763062]      + |              25     34765
##           tx_name
##       <character>
##   [1]  uc001uua.1
##   [2]  uc001uub.1
##   [3]  uc004bzv.3
##   [4]  uc004bzw.3
##   ---
##   seqlengths:
##                    chr1                 chr2 ...       chrUn_gl000249
##               249250621            243199373 ...                38502
```


The <emph>Rsamtools</emph> vignette illustrates a procedure for
creating a focused transcript set based on EBI biomart queries.
Here we obtain Entrez IDs for genes related to caffeine metabolism,
translate the feature names to ENSEMBL, and then create the
transcript database.  This requires internet connectivity and
nontrivial data interchange,
so is suppressed.

```r
library(KEGG.db)
kid <- revmap(KEGGPATHID2NAME)[["Caffeine metabolism"]]
egid <- KEGGPATHID2EXTID[[sprintf("hsa%s", kid)]]
library(biomaRt)
mart <- useMart("ensembl", "hsapiens_gene_ensembl")
ensid <- getBM(c("ensembl_transcript_id"), filters = "entrezgene", values = egid, 
    mart = mart)
library(GenomicFeatures)
txdb.caf <- makeTranscriptDbFromBiomart(transcript_ids = ensid)
transcripts(txdb.caf)
```

Partial result:
<pre>
> transcripts(txdb.caf)
GRanges with 22 ranges and 2 metadata columns:
       seqnames               ranges strand   |     tx_id         tx_name
          <Rle>            <IRanges>  <Rle>   | <integer>     <character>
   [1]       15 [75041185, 75048543]      +   |        18 ENST00000343932
   [2]       19 [41594368, 41602099]      +   |        19 ENST00000330436
   [3]       19 [41349444, 41356352]      -   |        20 ENST00000301141
   [4]       19 [41381344, 41388657]      -   |        21 ENST00000291764
   [5]       19 [41381344, 41388657]      -   |        22 ENST00000301146
</pre>

Exercise: Generalize the code above so that transcript sets
for arbitrary user-selected KEGG pathways are returned.
Introduce an option to use KEGGREST to retrieve the pathway
constituents.

<h3>Gene sets</h3>

Any vector of identifiers can represent a gene set.  A more formal
approach to managing information on gene sets is given in the 
<emph>GSEABase</emph> package, which is intended as infrastructure
for gene set enrichment testing.

In the following code, we use built-in data extracted from Broad Institute
MSIGDB to obtain a ``set collection'' representing two cytobands.

```r
library(GSEABase)
```

```
## Loading required package: annotate Loading required package: graph
```

```r
fl <- system.file("extdata", "Broad.xml", package = "GSEABase")
gs <- getBroadSets(fl)
gs
```

```
## GeneSetCollection
##   names: chr5q23, chr16q24 (2 total)
##   unique identifiers: ZNF474, CCDC100, ..., TRAPPC2L (215 total)
##   types in collection:
##     geneIdType: SymbolIdentifier (1 total)
##     collectionType: BroadCollection (1 total)
```

```r
names(gs)
```

```
## [1] "chr5q23"  "chr16q24"
```

```r
sapply(gs, function(x) length(geneIds(x)))
```

```
## [1]  86 129
```


Let's obtain the transcript GRanges for all the genes in \Sexpr{names(gs)[1]}.
First we'll get the Entrez IDs.

```r
egs = gs[[1]]
geneIdType(egs)
```

```
## geneIdType: Symbol
```

```r
SymbolIdentifier("ENTREZID")
```

```
## geneIdType: Symbol (ENTREZID)
```

```r
geneIdType(egs) = EntrezIdentifier("org.Hs.eg.db")
egs
```

```
## setName: chr5q23 
## geneIds: 133923, 9542, ..., 728586 (total: 58)
## geneIdType: EntrezId (org.Hs.eg.db)
## collectionType: Broad
##   bcCategory: c1 (Positional)
##   bcSubCategory:  NA
## details: use 'details(object)'
```

Notice the loss of tokens after the transformation.
With these, we call into transcripts again:

```r
tsel2 = transcripts(txdb, vals = list(gene_id = unique(geneIds(egs))), col = c("gene_id", 
    "tx_id", "tx_name"))
length(tsel2)
```

```
## [1] 128
```

```r
tsel2[1:3]
```

```
## GRanges with 3 ranges and 3 metadata columns:
##       seqnames               ranges strand |         gene_id     tx_id
##          <Rle>            <IRanges>  <Rle> | <CharacterList> <integer>
##   [1]     chr5 [52285156, 52390609]      + |            3673     20026
##   [2]     chr5 [52285156, 52390609]      + |            3673     20027
##   [3]     chr5 [52285156, 52390609]      + |            3673     20028
##           tx_name
##       <character>
##   [1]  uc003joy.3
##   [2]  uc011cqa.2
##   [3]  uc011cqb.2
##   ---
##   seqlengths:
##                    chr1                 chr2 ...       chrUn_gl000249
##               249250621            243199373 ...                38502
```


Exercise.  Create a function that takes a GeneSet instance and returns the
ranges of all transcripts in the set.  Allow an option to reduce the
ranges to ``gene regions''.

Redistribution of Broad MSIGDB is not permitted, but you can register to get
your own copy.

\includegraphics{msigFly}

Once you've acquired, e.g., the XML for the MSIGDB 3.1, the <code>getBroadSets</code> function can be used to obtain
a formal GeneSetCollection object representing the whole collection.  The resulting object can then be saved
for reuse.

<pre>
bash-3.2$ ls -tl ~/Downloads/m*xml
-rw-r--r--@ 1 stvjc  staff  50648647 Mar 16 07:43 /Users/stvjc/Downloads/msigdb_v3.1.xml
# in R: is not fast
msig3.1 = getBroadSets("/Users/stvjc/Downloads/msigdb_v3.1.xml")
save(msig3.1, file="msig3.1.rda")
## check size
bash-3.2$ ls -tl msig3.1.rda
-rw-r--r--  1 stvjc  staff  6227937 Mar 16 07:49 msig3.1.rda
</pre>

Here are some basic exploratory tasks, which you can't carry out unless you've
registered and acquired the MSIGDB as noted above.

<pre>
> msig3.1
GeneSetCollection
  names: NUCLEOPLASM, EXTRINSIC_TO_PLASMA_MEMBRANE, ..., LEF1_UP.V1_UP (8513 total)
  unique identifiers: HNRNPK, XRCC6, ..., LOC650488 (31847 total)
  types in collection:
    geneIdType: SymbolIdentifier (1 total)
    collectionType: BroadCollection (1 total)
> names(msig3.1)[1:10]
 [1] "NUCLEOPLASM"                                
 [2] "EXTRINSIC_TO_PLASMA_MEMBRANE"               
 [3] "ORGANELLE_PART"                             
 [4] "CELL_PROJECTION_PART"                       
 [5] "CYTOPLASMIC_VESICLE_MEMBRANE"               
 [6] "GOLGI_MEMBRANE"                             
 [7] "MITOCHONDRIAL_OUTER_MEMBRANE"               
 [8] "ORGANELLAR_RIBOSOME"                        
 [9] "ORGANELLAR_SMALL_RIBOSOMAL_SUBUNIT"         
[10] "INTRINSIC_TO_ENDOPLASMIC_RETICULUM_MEMBRANE"
> details(msig3.1[[1]])
setName: NUCLEOPLASM 
geneIds: HNRNPK, XRCC6, ..., PRM2 (total: 279)
geneIdType: Symbol
collectionType: Broad
  bcCategory: c5 (GO)
  bcSubCategory:  CC
setIdentifier: M12345
description: Genes annotated by the GO term GO:0005654. That part 
    of the nuclear content other than the chromosomes or the nucleolus.
organism: Homo sapiens
pubMedIds: 
urls: file://Users/stvjc/Downloads/msigdb_v3.1.xml
      http://amigo.geneontology.org/cgi-bin/amigo/go.cgi?view=
           details&search_constraint=terms&depth=0&query=GO:0005654
contributor: Gene Ontology
setVersion: 0.0.1
creationDate: Sat Mar 16 07:45:23 2013
</pre>

<h3>Pathway and network models</h3>

By ``model'' here we mean the data structure that represents a
gene set and relations among genes.  Algebraic graphs (node and
edge sets) are a basic tool for representation of networks, and Bioconductor's
<emph>graph</emph> package endeavors to support fairly rich computations
on graphs.  The Boost Graph Library of graph algorithms (www.boost.org)
is interfaced via package <emph>RBGL</emph>.

The Kyoto Encyclopedia of Genes and Genomes is a familiar repository of
gene network models.  Some changes to its open status have occurred in the
past few years, and creation of KEGG.db has been complicated by redistribution
restrictions.  Note the message after <code>library(KEGG.db)</code>.

However, other forms of access may not be so hampered.  The REST API
is employed by <emph>KEGGREST</emph> (which will be available in the release coming
in April 2013).  The <emph>KEGGgraph</emph> package manipulates XML models
of KEGG pathways, and the two vignettes of that package should be studied
with care.  Here are some illustrations.


```r
library(KEGGgraph)
```

```
## Loading required package: XML
## 
## Attaching package: 'XML'
## 
## The following object is masked from 'package:graph':
## 
## addNode
```

```r
sfile <- system.file("extdata/hsa04010.xml", package = "KEGGgraph")
mapksp <- parseKGML(sfile)
mapksp
```

```
## KEGG Pathway
## [ Title ]: MAPK signaling pathway
## [ Name ]: path:hsa04010
## [ Organism ]: hsa
## [ Number ] :04010
## [ Image ] :http://www.genome.jp/kegg/pathway/hsa/hsa04010.gif
## [ Link ] :http://www.genome.jp/dbget-bin/show_pathway?hsa04010
## ------------------------------------------------------------
## Statistics:
## 	136 node(s)
## 	171 edge(s)
## 	0 reaction(s)
## ------------------------------------------------------------
```

```r
class(mapksp)
```

```
## [1] "KEGGPathway"
## attr(,"package")
## [1] "KEGGgraph"
```

```r
nodes(mapksp)[1:2]
```

```
## $`1`
## KEGG Node (Entry '1'):
## ------------------------------------------------------------
## [ displayName ]: C00338
## [ Name ]: cpd:C00338
## [ Type ]: compound
## [ Link ]: http://www.genome.jp/dbget-bin/www_bget?compound+C00338
## ------------------------------------------------------------
## 
## $`2`
## KEGG Node (Entry '2'):
## ------------------------------------------------------------
## [ displayName ]: RASGRF1, GRF1...
## [ Name ]: hsa:5923,hsa:5924
## [ Type ]: gene
## [ Link ]: http://www.genome.jp/dbget-bin/www_bget?hsa+5923+5924
## ------------------------------------------------------------
```

```r
edges(mapksp)[1:2]
```

```
## $relation
##   KEGG Edge (Type: PPrel):
## ------------------------------------------------------------
## [ Entry 1 ID ]: 47
## [ Entry 2 ID ]: 40
## [ Subtype ]: 
##   [ Subtype name ]: activation
##   [ Subtype value ]: -->
## ------------------------------------------------------------
## 
## $relation
##   KEGG Edge (Type: PPrel):
## ------------------------------------------------------------
## [ Entry 1 ID ]: 46
## [ Entry 2 ID ]: 40
## [ Subtype ]: 
##   [ Subtype name ]: activation
##   [ Subtype value ]: -->
## ------------------------------------------------------------
```


The pathway object is not a graph instance <emph>per se</emph>.  We continue:

```r
mapkg = KEGGpathway2Graph(mapksp)
mapkg
```

```
## A graphNEL graph with directed edges
## Number of Nodes = 265 
## Number of Edges = 876
```

```r
nodes(mapkg)[1:2]
```

```
## [1] "hsa:5923" "hsa:5924"
```


We can get a rough sense of the topology of the
pathway with Rgraphviz.

\includegraphics{mapk}

Exercise: Show that the gene with the largest out degree value is MAPK1.

Other pathway catalogs with Bioconductor interfaces are Reactome
(package <emph>reactome.db</emph>), and cMAP.

<h3>SNPs</h3>

Annotation from NCBI dbSNP is acquired approximately every six months and
stored in the SNPlocs.Hsapiens.dbSNP.[yyyymmdd] package.

I am not requiring that you download this package as it is very
large.  If you do it, the following is illustrative.


```r
library(SNPlocs.Hsapiens.dbSNP.20120608)
```


```r
s22 = getSNPlocs("ch22", as.GRanges = TRUE)
length(s22)
```


```r
s22[1:3]
```


It is also possible to query the VCF distributed by
NCBI using code like the following.  This allows targeted
retrieval using a GRanges.


```r
dbsvcf = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz"
```




```r
getSNPlocs2
```

```
## function(ch, start=NA, end=NA, vcfpath=dbsvcf) {
## #
## # use VCF as backend for getSNPlocs request; ch can be "ch[n]" as previously
## # with ch="ch1", start=400000, end=600000 you get 491 ranges in a few sec over web
## #
##  require(Homo.sapiens)
##  cl = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
##  ch = gsub("ch", "", ch)
##  ch = paste0("chr", ch)
##  if (!ch %in% names(cl)) stop(paste("can't get chrlength of ", ch, "in TxDb"))
##  ln = cl[ch]
##  sn = gsub("chr", "", ch)
##  if (is.na(start)) start = 1
##  if (is.na(end)) end=ln
##  p = ScanVcfParam(which=GRanges(seqnames=sn, IRanges(start, end)))
##  ans = scanVcf(dbsvcf, param=p)
##  gr = ans[[1]]$rowData
##  gr$REF = ans[[1]]$REF
##  gr$ALT = ans[[1]]$ALT
##  gr
## }
```



<h3>Homology databases</h3>

The inparanoid ortholog databases are the primary
tool for homology considerations with Bioconductor.
Because the references have not been updated since 2009,
it is an open question whether these should be regarded as definitive.
This code is not evaluated but is provided for illustration of the
potential usage.


```r
library(hom.Hs.inp.db)
library(hom.Mm.inp.db)
library(org.Mm.eg.db)
```


The following function uses old patterns of symbol retrieval
to map gene symbols between human and mouse.

```r
HMSymMap = function(hsym) {
    eg = get(hsym, revmap(org.Hs.egSYMBOL))
    ep = get(eg, org.Hs.egENSEMBLPROT)
    mps = na.omit(unlist(mget(ep, hom.Hs.inpMUSMU, ifnotfound = NA)))
    megs = na.omit(unlist(mget(mps, revmap(org.Mm.egENSEMBLPROT), ifnotfound = NA)))
    unlist(mget(megs, org.Mm.egSYMBOL))
}
```


Run the function to identify mouse genes homologous to SERPINA3.


```r
HMSymMap("SERPINA3")
```


<h3>Other annotation concepts</h3>

The <emph>AnnotationHub</emph>
is a collection lightly curated R serializations
of annotation resources.  Its use requires
internet connectivity; as of Sept. 16 2013
the following computations occur.

<pre>
> library(AnnotationHub)
> ah = AnnotationHub()
> ah
class: AnnotationHub 
length: 8674 
filters: none 
hubUrl: http://annotationhub.bioconductor.org/ah 
snapshotVersion: 2.13; snapshotDate: 2013-06-29
hubCache: /Users/stvjc/.AnnotationHub 
> ahn = names(ah)
> sahn = strsplit(ahn, "\\.")
> table(sapply(sahn, "[", 1))

     dbSNP    ensembl goldenpath 
       302       1673       6699 
> cumsum(.Last.value)
     dbSNP    ensembl goldenpath 
       302       1975       8674 
> ahn[1]
[1] "dbSNP.organisms.human_9606.VCF.ByChromosome.01.12156.ASW.RData"
> ahn[303]
[1] "ensembl.release.69.fasta.ailuropoda_melanoleuca.cdna.Ailuropoda_melanoleuca.ailMel1.69.cdna.all.fa.rz"
> ahn[1976]
[1] "goldenpath.ailMel1.database.blastHg18KG_0.0.1.RData"
> grep("Hmm", ahn, value=TRUE, ignore.case=TRUE)
 [1] "goldenpath.hg18.database.wgEncodeBroadHmm_0.0.1.RData"
 [2] "goldenpath.hg19.encodeDCC.wgEncodeRikenCage.wgEncodeRikenCageA549CellPapTssHmm.bedRnaElements_0.0.1.RData"            
 [3] "goldenpath.hg19.encodeDCC.wgEncodeRikenCage.wgEncodeRikenCageA549CytosolPapTssHmm.bedRnaElements_0.0.1.RData"         
 # ...

> ah[["dbSNP.organisms.human_9606.VCF.ByChromosome.01.12156.ASW.RData"]]
Loading required package: VariantAnnotation
Loading required package: GenomicRanges
Loading required package: XVector
Loading required package: Rsamtools
Loading required package: Biostrings

Attaching package: 'VariantAnnotation'

The following object is masked from 'package:base':

    tabulate

class: CollapsedVCF 
dim: 124907 49 
rowData(vcf):
  GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
info(vcf):
  DataFrame with 52 columns: AF, NS, RSPOS, RV, VP, GENEINFO, dbSNPBuildID, ...
info(header(vcf)):
                Number Type   
   AF           .      Float  
   NS           1      Integer
   RSPOS        1      Integer
   RV           0      Flag   
   VP           1      String 
   GENEINFO     1      String 
   dbSNPBuildID 1      Integer
   SAO          1      Integer
   SSR          1      Integer
   GMAF         1      Float  
   WGT          1      Integer
   VC           1      String 
   PM           0      Flag   
   TPA          0      Flag   
   PMC          0      Flag   
   S3D          0      Flag   
   SLO          0      Flag   
   NSF          0      Flag   
   NSM          0      Flag   
   NSN          0      Flag   
   REF          0      Flag   
   SYN          0      Flag   
   U3           0      Flag   
   U5           0      Flag   
   ASS          0      Flag   
   DSS          0      Flag   
   INT          0      Flag   
   R3           0      Flag   
   R5           0      Flag   
   OTH          0      Flag   
   CFL          0      Flag   
   ASP          0      Flag   
   MUT          0      Flag   
   VLD          0      Flag   
   G5A          0      Flag   
   G5           0      Flag   
   HD           0      Flag   
   GNO          0      Flag   
   KGValidated  0      Flag   
   KGPhase1     0      Flag   
   KGPilot123   0      Flag   
   KGPROD       0      Flag   
   OTHERKG      0      Flag   
   PH3          0      Flag   
   CDA          0      Flag   
   LSD          0      Flag   
   MTP          0      Flag   
   OM           0      Flag   
   NOC          0      Flag   
   WTD          0      Flag   
   NOV          0      Flag   
   GCF          0      Flag   
                Description                                       
   AF           Allele Frequency, for each ALT allele, in the s...
   NS           Number of Samples With Data                       
   RSPOS        Chr position reported in dbSNP                    
   RV           RS orientation is reversed                        
   VP           Variation Property.  Documentation is at ftp://...   GENEINFO     Pairs each of gene symbol:gene id.  The gene sy...
   dbSNPBuildID First dbSNP Build for RS                          
   SAO          Variant Allele Origin: 0 - unspecified, 1 - Ger...
   SSR          Variant Suspect Reason Code, 0 - unspecified, 1...
   GMAF         Global Minor Allele Frequency [0, 0.5]; global ...
   WGT          Weight, 00 - unmapped, 1 - weight 1, 2 - weight...
   VC           Variation Class                                   
   PM           Variant is Precious(Clinical,Pubmed Cited)        
   TPA          Provisional Third Party Annotation(TPA) (curren...
   PMC          Links exist to PubMed Central article             
   S3D          Has 3D structure - SNP3D table                    
   SLO          Has SubmitterLinkOut - From SNP->SubSNP->Batch....
   NSF          Has non-synonymous frameshift A coding region v...
   NSM          Has non-synonymous missense A coding region var...
   NSN          Has non-synonymous nonsense A coding region var...
   REF          Has reference A coding region variation where o...
   SYN          Has synonymous A coding region variation where ...
   U3           In 3' UTR Location is in an untranslated region...
   U5           In 5' UTR Location is in an untranslated region...
   ASS          In acceptor splice site FxnCode = 73              
   DSS          In donor splice-site FxnCode = 75                 
   INT          In Intron FxnCode = 6                             
   R3           In 3' gene region FxnCode = 13                    
   R5           In 5' gene region FxnCode = 15                    
   OTH          Has other variant with exactly the same set of ...
   CFL          Has Assembly conflict. This is for weight 1 and...
   ASP          Is Assembly specific. This is set if the varian...
   MUT          Is mutation (journal citation, explicit fact): ...
   VLD          Is Validated.  This bit is set if the variant h...
   G5A          >5% minor allele frequency in each and all popu...
   G5           >5% minor allele frequency in 1+ populations      
   HD           Marker is on high density genotyping kit (50K d...
   GNO          Genotypes available. The variant has individual...
   KGValidated  1000 Genome validated                             
   KGPhase1     1000 Genome phase 1 (incl. June Interim phase 1)  
   KGPilot123   1000 Genome discovery all pilots 2010(1,2,3)      
   KGPROD       Has 1000 Genome submission                        
   OTHERKG      non-1000 Genome submission                        
   PH3          HAP_MAP Phase 3 genotyped: filtered, non-redundant
   CDA          Variation is interrogated in a clinical diagnos...
   LSD          Submitted from a locus-specific database          
   MTP          Microattribution/third-party annotation(TPA:GWA...
   OM           Has OMIM/OMIA                                     
   NOC          Contig allele not present in variant allele lis...
   WTD          Is Withdrawn by submitter If one member ss is w...
   NOV          Rs cluster has non-overlapping allele sets. Tru...
   GCF          Has Genotype Conflict Same (rs, ind), different...
geno(vcf):
  SimpleList of length 6: GT, GT2, GT3, GC, GC2, GC3
geno(header(vcf)):
       Number Type    Description          
   GT  1      String  Genotype             
   GT2 1      String  Second Genotype      
   GT3 1      String  Third Genotype       
   GC  1      Integer Genotype Count       
   GC2 1      Integer Second Genotype Count
   GC3 1      Integer Third Genotype Count 
</pre>

