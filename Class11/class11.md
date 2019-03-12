Class 11: Structural Bioinformatics
================
Paige Dubelko
February 14, 2019

Section 1. The PDB
------------------

``` r
pdb <- read.csv("Data Export Summary.csv", row.names = 1)
```

Let's determine the percentage of structures solved by X-Ray and Electron Microscopy.

``` r
xRayPercent <- (pdb$Total[1]/ sum(pdb$Total))*100
xRayPercent
```

    ## [1] 89.43069

``` r
emPercent <- (pdb$Total[3]/sum(pdb$Total))*100
emPercent
```

    ## [1] 1.88697

Another way to do this is:

``` r
precent.by.method <- (pdb$Total/sum(pdb$Total))*100
names(precent.by.method) <- rownames(pdb)
precent.by.method
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.43068692          8.40846082          1.88696977 
    ##               Other        Multi Method 
    ##          0.18325960          0.09062288

From the website, what proportion of structures are protein?

``` r
proteinPercent <- sum(pdb$Proteins)/sum(pdb$Total)*100
proteinPercent
```

    ## [1] 92.76561

Getting started with Bio3D
--------------------------

``` r
library(bio3d)
```

``` r
hsg <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
head(aa321(hsg$atom$resid))
```

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: MK1

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## Warning in FUN(X[[i]], ...): Unknown 3-letters code for aminoacid: HOH

    ## [1] "P" "P" "P" "P" "P" "P"

``` r
prot.inds <- atom.select(hsg, string = 'protein')
prot.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = hsg, string = "protein")
    ## 
    ##    Atom Indices#: 1514  ($atom)
    ##    XYZ  Indices#: 4542  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
lig.inds <- atom.select(hsg, string = 'ligand')
```

``` r
prot.pdb <- trim.pdb(hsg, prot.inds)
write.pdb(prot.pdb, file = "protein.pdb")
lig.pdb <- trim.pdb(hsg, lig.inds)
write.pdb(lig.pdb, file = "ligand.pdb")
```
