Class 18: Cancer Genomics
================
Paige Dubelko

``` r
library(bio3d)
```

Read in FASTA

``` r
alignment <- read.fasta("lecture18_sequences.fa")
alignment
```

    ##              1        .         .         .         .         .         60 
    ## P53_wt       MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    ## P53_mutant   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
    ##              **************************************** ******************* 
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## P53_wt       DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ## P53_mutant   DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ##              **** ******************************************************* 
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## P53_wt       SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ## P53_mutant   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ##              ************************************************************ 
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .         .         .         240 
    ## P53_wt       RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ## P53_mutant   RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ##              ******************************** *************************** 
    ##            181        .         .         .         .         .         240 
    ## 
    ##            241        .         .         .         .         .         300 
    ## P53_wt       SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    ## P53_mutant   SCMGGMNRRPILTIITLEV-----------------------------------------
    ##              ******************                                           
    ##            241        .         .         .         .         .         300 
    ## 
    ##            301        .         .         .         .         .         360 
    ## P53_wt       PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    ## P53_mutant   ------------------------------------------------------------
    ##                                                                           
    ##            301        .         .         .         .         .         360 
    ## 
    ##            361        .         .         .  393 
    ## P53_wt       GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    ## P53_mutant   ---------------------------------
    ##                                                
    ##            361        .         .         .  393 
    ## 
    ## Call:
    ##   read.fasta(file = "lecture18_sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 393 position columns (259 non-gap, 134 gap) 
    ## 
    ## + attr: id, ali, call

Use bio3d package to align fasta file

``` r
names <-  rownames(alignment)
alignedSeq <- seqaln(alignment)
```

Identify which positions have mismatches

``` r
ide <- conserv(alignedSeq$ali, method = "identity")
mutant.sites <- which(ide < 1)

gaps <- gap.inspect(alignedSeq)
mutant.sites <- mutant.sites[mutant.sites %in% gaps$f.inds]

mutant.sites
```

    ## [1]  41  65 213 259

Now extract all 9-mer mutants encompassing sequences for each mutant site. This is equivalent to finding the sequence region eight residues before and eight residues after our mutation sites.

``` r
## Make a "names" label for our output sequences (one per mutant)
mutant.names <- paste0(alignedSeq$ali["P53_wt",mutant.sites],
                       mutant.sites,
                       alignedSeq$ali["P53_mutant",mutant.sites])

mutant.names
```

    ## [1] "D41L"  "R65W"  "R213V" "D259V"

``` r
## Sequence positions surounding each mutant site
start.position <- mutant.sites - 8
end.position <-  mutant.sites + 8

# Blank matrix to store sub-sequences
store.seqs <- matrix("-", nrow=length(mutant.sites), ncol=17)
rownames(store.seqs) <- mutant.names

## Extract each sub-sequence
for(i in 1:length(mutant.sites)) {
  store.seqs[i,] <- alignedSeq$ali["P53_mutant",start.position[i]:end.position[i]]
}

store.seqs
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## D41L  "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"   "L"  
    ## R65W  "D"  "P"  "G"  "P"  "D"  "E"  "A"  "P"  "W"  "M"   "P"   "E"   "A"  
    ## R213V "Y"  "L"  "D"  "D"  "R"  "N"  "T"  "F"  "V"  "H"   "S"   "V"   "V"  
    ## D259V "I"  "L"  "T"  "I"  "I"  "T"  "L"  "E"  "V"  "-"   "-"   "-"   "-"  
    ##       [,14] [,15] [,16] [,17]
    ## D41L  "S"   "P"   "D"   "D"  
    ## R65W  "A"   "P"   "P"   "V"  
    ## R213V "V"   "P"   "Y"   "E"  
    ## D259V "-"   "-"   "-"   "-"

``` r
#output FASTA file
write.fasta(seqs = store.seqs, id = mutant.names, file = "subsequences.fa")
```

``` r
result <- read.csv("result.csv")
result
```

    ##          allele seq_num start end length   peptide
    ## 1   HLA-A*02:01       3     1   9      9 YLDDRNTFV
    ## 2   HLA-B*35:01       3     8  16      9 FVHSVVVPY
    ## 3   HLA-B*07:02       1     1   9      9 SPLPSQAML
    ## 4   HLA-A*02:01       2     9  17      9 WMPEAAPPV
    ## 5   HLA-B*07:02       1     3  11      9 LPSQAMLDL
    ## 6   HLA-A*02:01       4     1   9      9 ILTIITLEV
    ## 7   HLA-A*68:01       3     8  16      9 FVHSVVVPY
    ## 8   HLA-B*35:01       1     3  11      9 LPSQAMLDL
    ## 9   HLA-B*35:01       2     7  15      9 APWMPEAAP
    ## 10  HLA-B*07:02       2     7  15      9 APWMPEAAP
    ## 11  HLA-A*02:01       1     7  15      9 AMLDLMLSP
    ## 12  HLA-B*35:01       2     6  14      9 EAPWMPEAA
    ## 13  HLA-B*35:01       1     1   9      9 SPLPSQAML
    ## 14  HLA-B*35:01       2     1   9      9 DPGPDEAPW
    ## 15  HLA-B*35:01       2     3  11      9 GPDEAPWMP
    ## 16  HLA-A*02:01       1     5  13      9 SQAMLDLML
    ## 17  HLA-A*02:01       3     6  14      9 NTFVHSVVV
    ## 18  HLA-A*02:01       1     8  16      9 MLDLMLSPD
    ## 19  HLA-B*35:01       1     6  14      9 QAMLDLMLS
    ## 20  HLA-B*35:01       1     8  16      9 MLDLMLSPD
    ## 21  HLA-A*68:01       3     6  14      9 NTFVHSVVV
    ## 22  HLA-B*07:02       3     5  13      9 RNTFVHSVV
    ## 23  HLA-B*35:01       2     9  17      9 WMPEAAPPV
    ## 24  HLA-B*07:02       1     5  13      9 SQAMLDLML
    ## 25  HLA-B*07:02       2     1   9      9 DPGPDEAPW
    ## 26  HLA-A*02:01       3     8  16      9 FVHSVVVPY
    ## 27  HLA-B*35:01       3     2  10      9 LDDRNTFVH
    ## 28  HLA-B*35:01       1     4  12      9 PSQAMLDLM
    ## 29  HLA-B*07:02       3     1   9      9 YLDDRNTFV
    ## 30  HLA-B*07:02       4     1   9      9 ILTIITLEV
    ## 31  HLA-B*07:02       3     6  14      9 NTFVHSVVV
    ## 32  HLA-A*68:01       1     6  14      9 QAMLDLMLS
    ## 33  HLA-B*07:02       3     8  16      9 FVHSVVVPY
    ## 34  HLA-B*07:02       2     9  17      9 WMPEAAPPV
    ## 35  HLA-B*35:01       2     8  16      9 PWMPEAAPP
    ## 36  HLA-A*68:01       2     6  14      9 EAPWMPEAA
    ## 37  HLA-B*35:01       3     6  14      9 NTFVHSVVV
    ## 38  HLA-A*02:01       1     6  14      9 QAMLDLMLS
    ## 39  HLA-A*68:01       3     7  15      9 TFVHSVVVP
    ## 40  HLA-B*35:01       1     5  13      9 SQAMLDLML
    ## 41  HLA-B*35:01       1     9  17      9 LDLMLSPDD
    ## 42  HLA-B*07:02       1     7  15      9 AMLDLMLSP
    ## 43  HLA-A*68:01       3     9  17      9 VHSVVVPYE
    ## 44  HLA-B*35:01       3     7  15      9 TFVHSVVVP
    ## 45  HLA-A*02:01       2     5  13      9 DEAPWMPEA
    ## 46  HLA-A*68:01       4     1   9      9 ILTIITLEV
    ## 47  HLA-A*02:01       1     1   9      9 SPLPSQAML
    ## 48  HLA-B*07:02       2     3  11      9 GPDEAPWMP
    ## 49  HLA-B*07:02       2     6  14      9 EAPWMPEAA
    ## 50  HLA-A*02:01       3     4  12      9 DRNTFVHSV
    ## 51  HLA-B*35:01       2     2  10      9 PGPDEAPWM
    ## 52  HLA-B*35:01       4     1   9      9 ILTIITLEV
    ## 53  HLA-B*07:02       1     8  16      9 MLDLMLSPD
    ## 54  HLA-A*02:01       1     3  11      9 LPSQAMLDL
    ## 55  HLA-B*07:02       1     4  12      9 PSQAMLDLM
    ## 56  HLA-A*02:01       3     5  13      9 RNTFVHSVV
    ## 57  HLA-A*68:01       3     1   9      9 YLDDRNTFV
    ## 58  HLA-B*35:01       1     7  15      9 AMLDLMLSP
    ## 59  HLA-A*68:01       3     5  13      9 RNTFVHSVV
    ## 60  HLA-A*02:01       1     2  10      9 PLPSQAMLD
    ## 61  HLA-A*02:01       2     3  11      9 GPDEAPWMP
    ## 62  HLA-B*35:01       2     5  13      9 DEAPWMPEA
    ## 63  HLA-B*07:02       3     2  10      9 LDDRNTFVH
    ## 64  HLA-B*07:02       2     8  16      9 PWMPEAAPP
    ## 65  HLA-A*02:01       1     4  12      9 PSQAMLDLM
    ## 66  HLA-A*68:01       1     8  16      9 MLDLMLSPD
    ## 67  HLA-A*02:01       2     6  14      9 EAPWMPEAA
    ## 68  HLA-A*68:01       2     5  13      9 DEAPWMPEA
    ## 69  HLA-B*35:01       3     4  12      9 DRNTFVHSV
    ## 70  HLA-B*07:02       1     6  14      9 QAMLDLMLS
    ## 71  HLA-B*35:01       3     1   9      9 YLDDRNTFV
    ## 72  HLA-A*68:01       1     5  13      9 SQAMLDLML
    ## 73  HLA-A*68:01       3     2  10      9 LDDRNTFVH
    ## 74  HLA-B*07:02       2     4  12      9 PDEAPWMPE
    ## 75  HLA-A*02:01       3     9  17      9 VHSVVVPYE
    ## 76  HLA-A*02:01       2     8  16      9 PWMPEAAPP
    ## 77  HLA-A*68:01       1     1   9      9 SPLPSQAML
    ## 78  HLA-A*02:01       3     7  15      9 TFVHSVVVP
    ## 79  HLA-A*02:01       2     2  10      9 PGPDEAPWM
    ## 80  HLA-A*02:01       2     7  15      9 APWMPEAAP
    ## 81  HLA-A*68:01       1     4  12      9 PSQAMLDLM
    ## 82  HLA-B*35:01       2     4  12      9 PDEAPWMPE
    ## 83  HLA-A*68:01       1     7  15      9 AMLDLMLSP
    ## 84  HLA-A*02:01       1     9  17      9 LDLMLSPDD
    ## 85  HLA-A*68:01       2     9  17      9 WMPEAAPPV
    ## 86  HLA-B*07:02       3     4  12      9 DRNTFVHSV
    ## 87  HLA-A*02:01       3     2  10      9 LDDRNTFVH
    ## 88  HLA-A*68:01       3     4  12      9 DRNTFVHSV
    ## 89  HLA-B*07:02       1     9  17      9 LDLMLSPDD
    ## 90  HLA-A*68:01       3     3  11      9 DDRNTFVHS
    ## 91  HLA-B*07:02       2     5  13      9 DEAPWMPEA
    ## 92  HLA-B*35:01       3     5  13      9 RNTFVHSVV
    ## 93  HLA-A*68:01       2     1   9      9 DPGPDEAPW
    ## 94  HLA-B*07:02       2     2  10      9 PGPDEAPWM
    ## 95  HLA-B*35:01       3     3  11      9 DDRNTFVHS
    ## 96  HLA-B*07:02       3     3  11      9 DDRNTFVHS
    ## 97  HLA-A*68:01       1     3  11      9 LPSQAMLDL
    ## 98  HLA-A*68:01       2     4  12      9 PDEAPWMPE
    ## 99  HLA-A*68:01       2     7  15      9 APWMPEAAP
    ## 100 HLA-B*35:01       3     9  17      9 VHSVVVPYE
    ## 101 HLA-A*68:01       1     2  10      9 PLPSQAMLD
    ## 102 HLA-A*68:01       2     8  16      9 PWMPEAAPP
    ## 103 HLA-A*68:01       2     3  11      9 GPDEAPWMP
    ## 104 HLA-B*07:02       3     7  15      9 TFVHSVVVP
    ## 105 HLA-A*02:01       2     1   9      9 DPGPDEAPW
    ## 106 HLA-A*68:01       1     9  17      9 LDLMLSPDD
    ## 107 HLA-B*07:02       3     9  17      9 VHSVVVPYE
    ## 108 HLA-A*02:01       2     4  12      9 PDEAPWMPE
    ## 109 HLA-A*02:01       3     3  11      9 DDRNTFVHS
    ## 110 HLA-B*07:02       1     2  10      9 PLPSQAMLD
    ## 111 HLA-B*35:01       1     2  10      9 PLPSQAMLD
    ## 112 HLA-A*68:01       2     2  10      9 PGPDEAPWM
    ##                                     method percentile_rank ann_ic50
    ## 1   Consensus (ann/comblib_sidney2008/smm)             0.2     4.05
    ## 2   Consensus (ann/comblib_sidney2008/smm)             0.2     5.41
    ## 3   Consensus (ann/comblib_sidney2008/smm)             0.4    43.33
    ## 4   Consensus (ann/comblib_sidney2008/smm)             0.4     5.78
    ## 5   Consensus (ann/comblib_sidney2008/smm)             0.5    35.11
    ## 6   Consensus (ann/comblib_sidney2008/smm)             0.7    35.81
    ## 7                      Consensus (ann/smm)             2.1    74.54
    ## 8   Consensus (ann/comblib_sidney2008/smm)             2.2    67.78
    ## 9   Consensus (ann/comblib_sidney2008/smm)             2.2  1940.03
    ## 10  Consensus (ann/comblib_sidney2008/smm)             2.8  1648.70
    ## 11  Consensus (ann/comblib_sidney2008/smm)             4.0   642.57
    ## 12  Consensus (ann/comblib_sidney2008/smm)             4.2  5013.09
    ## 13  Consensus (ann/comblib_sidney2008/smm)             5.2  1071.67
    ## 14  Consensus (ann/comblib_sidney2008/smm)             5.2  1736.36
    ## 15  Consensus (ann/comblib_sidney2008/smm)             6.6  8634.42
    ## 16  Consensus (ann/comblib_sidney2008/smm)             7.7  1075.00
    ## 17  Consensus (ann/comblib_sidney2008/smm)             8.5  2005.42
    ## 18  Consensus (ann/comblib_sidney2008/smm)             8.9  8067.39
    ## 19  Consensus (ann/comblib_sidney2008/smm)             9.2 12339.15
    ## 20  Consensus (ann/comblib_sidney2008/smm)             9.8 13127.78
    ## 21                     Consensus (ann/smm)             9.9  4337.95
    ## 22  Consensus (ann/comblib_sidney2008/smm)            12.0 12125.27
    ## 23  Consensus (ann/comblib_sidney2008/smm)            13.0  6749.78
    ## 24  Consensus (ann/comblib_sidney2008/smm)            16.0 15900.99
    ## 25  Consensus (ann/comblib_sidney2008/smm)            16.0 18475.22
    ## 26  Consensus (ann/comblib_sidney2008/smm)            16.0  7431.44
    ## 27  Consensus (ann/comblib_sidney2008/smm)            17.0 11962.78
    ## 28  Consensus (ann/comblib_sidney2008/smm)            18.0 17855.96
    ## 29  Consensus (ann/comblib_sidney2008/smm)            20.0 21242.54
    ## 30  Consensus (ann/comblib_sidney2008/smm)            20.0 18572.62
    ## 31  Consensus (ann/comblib_sidney2008/smm)            21.0 19835.71
    ## 32                     Consensus (ann/smm)            21.5 16986.06
    ## 33  Consensus (ann/comblib_sidney2008/smm)            22.0 16491.48
    ## 34  Consensus (ann/comblib_sidney2008/smm)            23.0 21201.89
    ## 35  Consensus (ann/comblib_sidney2008/smm)            23.0 23724.50
    ## 36                     Consensus (ann/smm)            24.5 22317.85
    ## 37  Consensus (ann/comblib_sidney2008/smm)            25.0 13148.96
    ## 38  Consensus (ann/comblib_sidney2008/smm)            26.0 12700.95
    ## 39                     Consensus (ann/smm)            27.0 24109.80
    ## 40  Consensus (ann/comblib_sidney2008/smm)            27.0 11363.37
    ## 41  Consensus (ann/comblib_sidney2008/smm)            28.0 24432.01
    ## 42  Consensus (ann/comblib_sidney2008/smm)            29.0 27382.15
    ## 43                     Consensus (ann/smm)            30.0 32555.75
    ## 44  Consensus (ann/comblib_sidney2008/smm)            31.0 24211.76
    ## 45  Consensus (ann/comblib_sidney2008/smm)            32.0 32574.78
    ## 46                     Consensus (ann/smm)            32.5 21351.99
    ## 47  Consensus (ann/comblib_sidney2008/smm)            35.0 22832.49
    ## 48  Consensus (ann/comblib_sidney2008/smm)            36.0 30434.90
    ## 49  Consensus (ann/comblib_sidney2008/smm)            36.0 30222.26
    ## 50  Consensus (ann/comblib_sidney2008/smm)            36.0 27305.53
    ## 51  Consensus (ann/comblib_sidney2008/smm)            36.0 22780.17
    ## 52  Consensus (ann/comblib_sidney2008/smm)            36.0 20045.64
    ## 53  Consensus (ann/comblib_sidney2008/smm)            41.0 23760.20
    ## 54  Consensus (ann/comblib_sidney2008/smm)            42.0 23975.82
    ## 55  Consensus (ann/comblib_sidney2008/smm)            43.0 25324.78
    ## 56  Consensus (ann/comblib_sidney2008/smm)            43.0 27402.61
    ## 57                     Consensus (ann/smm)            45.0 30003.32
    ## 58  Consensus (ann/comblib_sidney2008/smm)            45.0 27148.80
    ## 59                     Consensus (ann/smm)            46.0 34078.52
    ## 60  Consensus (ann/comblib_sidney2008/smm)            47.0 36359.52
    ## 61  Consensus (ann/comblib_sidney2008/smm)            47.0 33777.88
    ## 62  Consensus (ann/comblib_sidney2008/smm)            47.0 21641.09
    ## 63  Consensus (ann/comblib_sidney2008/smm)            48.0 30762.99
    ## 64  Consensus (ann/comblib_sidney2008/smm)            49.0 34187.83
    ## 65  Consensus (ann/comblib_sidney2008/smm)            49.0 27523.54
    ## 66                     Consensus (ann/smm)            50.5 34201.90
    ## 67  Consensus (ann/comblib_sidney2008/smm)            51.0 30525.92
    ## 68                     Consensus (ann/smm)            51.0 38144.75
    ## 69  Consensus (ann/comblib_sidney2008/smm)            53.0 35430.21
    ## 70  Consensus (ann/comblib_sidney2008/smm)            54.0 23185.75
    ## 71  Consensus (ann/comblib_sidney2008/smm)            55.0 25722.72
    ## 72                     Consensus (ann/smm)            56.5 34573.96
    ## 73                     Consensus (ann/smm)            57.0 38533.01
    ## 74  Consensus (ann/comblib_sidney2008/smm)            58.0 34956.88
    ## 75  Consensus (ann/comblib_sidney2008/smm)            58.0 30933.88
    ## 76  Consensus (ann/comblib_sidney2008/smm)            59.0 41313.26
    ## 77                     Consensus (ann/smm)            60.0 38342.53
    ## 78  Consensus (ann/comblib_sidney2008/smm)            61.0 30820.96
    ## 79  Consensus (ann/comblib_sidney2008/smm)            63.0 40729.61
    ## 80  Consensus (ann/comblib_sidney2008/smm)            63.0 32672.19
    ## 81                     Consensus (ann/smm)            64.0 39086.86
    ## 82  Consensus (ann/comblib_sidney2008/smm)            64.0 38130.71
    ## 83                     Consensus (ann/smm)            64.5 37852.43
    ## 84  Consensus (ann/comblib_sidney2008/smm)            66.0 40560.75
    ## 85                     Consensus (ann/smm)            67.5 38265.04
    ## 86  Consensus (ann/comblib_sidney2008/smm)            69.0 39889.64
    ## 87  Consensus (ann/comblib_sidney2008/smm)            69.0 38658.71
    ## 88                     Consensus (ann/smm)            69.5 41509.07
    ## 89  Consensus (ann/comblib_sidney2008/smm)            70.0 38207.93
    ## 90                     Consensus (ann/smm)            70.0 41664.30
    ## 91  Consensus (ann/comblib_sidney2008/smm)            73.0 32692.00
    ## 92  Consensus (ann/comblib_sidney2008/smm)            74.0 31115.15
    ## 93                     Consensus (ann/smm)            74.5 43473.68
    ## 94  Consensus (ann/comblib_sidney2008/smm)            75.0 36700.62
    ## 95  Consensus (ann/comblib_sidney2008/smm)            75.0 40176.83
    ## 96  Consensus (ann/comblib_sidney2008/smm)            76.0 39201.22
    ## 97                     Consensus (ann/smm)            76.0 43401.77
    ## 98                     Consensus (ann/smm)            78.5 45876.93
    ## 99                     Consensus (ann/smm)            79.0 44743.05
    ## 100 Consensus (ann/comblib_sidney2008/smm)            81.0 36486.03
    ## 101                    Consensus (ann/smm)            82.0 42850.70
    ## 102                    Consensus (ann/smm)            85.0 44392.03
    ## 103                    Consensus (ann/smm)            85.5 46076.43
    ## 104 Consensus (ann/comblib_sidney2008/smm)            86.0 35036.77
    ## 105 Consensus (ann/comblib_sidney2008/smm)            86.0 39659.42
    ## 106                    Consensus (ann/smm)            87.0 43683.03
    ## 107 Consensus (ann/comblib_sidney2008/smm)            88.0 31508.82
    ## 108 Consensus (ann/comblib_sidney2008/smm)            89.0 40690.42
    ## 109 Consensus (ann/comblib_sidney2008/smm)            90.0 37589.58
    ## 110 Consensus (ann/comblib_sidney2008/smm)            91.0 38680.45
    ## 111 Consensus (ann/comblib_sidney2008/smm)            92.0 42432.24
    ## 112                    Consensus (ann/smm)            98.0 46779.17
    ##     ann_rank    smm_ic50 smm_rank comblib_sidney2008_score
    ## 1       0.02        3.90      0.2                 3.63e-06
    ## 2       0.02        6.16      0.2                 0.000698
    ## 3       0.18       50.11      0.4                 0.000137
    ## 4       0.05       17.66      0.4                 9.93e-06
    ## 5       0.14       66.51      0.5                 6.56e-05
    ## 6       0.39       41.03      0.7                 1.91e-05
    ## 7       0.61      273.84      3.6                        -
    ## 8       0.24      170.44      2.2                   0.0026
    ## 9       2.20      710.53      5.8                 0.000187
    ## 10      2.80      899.27      3.1                 3.26e-05
    ## 11      3.50      358.19      4.0                 0.000736
    ## 12      4.20     1195.58      8.3                 9.91e-05
    ## 13      1.50      592.35      5.2                  0.00267
    ## 14      2.10      589.63      5.2                   0.0024
    ## 15      6.60     5100.12     21.0                 0.000598
    ## 16      4.70     1104.36      7.7                     0.01
    ## 17      6.90     1315.56      8.5                  0.00124
    ## 18     18.00     1429.26      8.9                 3.44e-06
    ## 19      9.20     7287.52     24.0                 0.000861
    ## 20      9.80     2521.04     14.0                 0.000523
    ## 21      6.80     1627.38     13.0                        -
    ## 22     10.00     6395.73     12.0                  0.00054
    ## 23      5.30     1956.95     13.0                     0.01
    ## 24     14.00    11063.44     16.0                  0.00219
    ## 25     16.00     6868.94     12.0                  0.00448
    ## 26     17.00     3707.75     16.0                 2.29e-05
    ## 27      9.00     3504.13     17.0                  0.00443
    ## 28     16.00    11955.84     31.0                  0.00199
    ## 29     19.00    18488.00     20.0                     0.01
    ## 30     16.00    19226.05     20.0                     0.01
    ## 31     18.00    20696.17     21.0                     0.01
    ## 32     18.00     6345.85     25.0                        -
    ## 33     15.00    22227.47     22.0                     0.01
    ## 34     19.00    23817.16     23.0                  0.00259
    ## 35     23.00    50070.28     52.0                  0.00273
    ## 36     23.00     7038.66     26.0                        -
    ## 37      9.90    14111.72     33.0                  0.00313
    ## 38     24.00    12885.76     26.0                 0.000853
    ## 39     25.00     8902.05     29.0                        -
    ## 40      8.50     9153.35     27.0                     0.03
    ## 41     24.00     9540.69     28.0                     0.01
    ## 42     29.00   206485.71     51.0                 0.000733
    ## 43     38.00     4925.95     22.0                        -
    ## 44     24.00    12290.80     31.0                     0.01
    ## 45     63.00    25888.69     32.0                 0.000311
    ## 46     22.00    23631.43     43.0                        -
    ## 47     40.00    13004.99     26.0                  0.00104
    ## 48     36.00    88489.15     39.0                  0.00197
    ## 49     36.00    92021.65     39.0                 0.000418
    ## 50     49.00    35327.26     36.0                 0.000742
    ## 51     22.00    18263.33     36.0                     0.03
    ## 52     18.00    21507.01     39.0                     0.01
    ## 53     23.00   322767.65     57.0                  0.00252
    ## 54     42.00    44989.38     39.0                     0.01
    ## 55     25.00   118274.19     43.0                     0.02
    ## 56     49.00    55476.62     41.0                  0.00148
    ## 57     33.00    52178.31     57.0                        -
    ## 58     29.00    32777.82     45.0                     0.03
    ## 59     42.00    36015.55     50.0                        -
    ## 60     75.00    92920.17     47.0                 0.000262
    ## 61     66.00    93134.37     47.0                  0.00157
    ## 62     20.00    35692.75     47.0                     0.03
    ## 63     37.00   174940.35     48.0                     0.02
    ## 64     49.00  2511250.29     84.0                  0.00211
    ## 65     49.00   219336.04     58.0                 0.000553
    ## 66     42.00    61587.13     59.0                        -
    ## 67     57.00   123626.05     51.0                  0.00111
    ## 68     53.00    32247.05     49.0                        -
    ## 69     53.00   263985.46     80.0                  0.00388
    ## 70     22.00   252284.17     54.0                     0.01
    ## 71     27.00    56439.09     55.0                     0.04
    ## 72     43.00   117622.39     70.0                        -
    ## 73     55.00    59770.93     59.0                        -
    ## 74     52.00  1819240.02     81.0                     0.01
    ## 75     58.00   796361.03     74.0                 0.000773
    ## 76     91.00   232332.52     59.0                 9.36e-05
    ## 77     54.00    90258.89     66.0                        -
    ## 78     57.00   285831.44     61.0                     0.01
    ## 79     90.00   347624.20     63.0                 0.000546
    ## 80     63.00   420833.21     66.0                  0.00298
    ## 81     57.00   119258.68     71.0                        -
    ## 82     64.00   244669.59     79.0                  0.00212
    ## 83     52.00   181338.49     77.0                        -
    ## 84     89.00   420833.21     66.0                  0.00325
    ## 85     54.00   231467.49     81.0                        -
    ## 86     81.00   719266.78     69.0                 0.000623
    ## 87     83.00   549680.08     69.0                 0.000961
    ## 88     67.00   128970.39     72.0                        -
    ## 89     70.00 21227067.45     99.0                     0.01
    ## 90     68.00   130463.79     72.0                        -
    ## 91     43.00  1261507.97     76.0                     0.01
    ## 92     39.00   182633.31     74.0                     0.04
    ## 93     78.00   124305.40     71.0                        -
    ## 94     61.00  1153161.14     75.0                     0.04
    ## 95     75.00   385993.51     85.0                     0.01
    ## 96     76.00  2454087.26     84.0                     0.01
    ## 97     77.00   161618.10     75.0                        -
    ## 98     92.00    83270.28     65.0                        -
    ## 99     85.00   134428.08     73.0                        -
    ## 100    57.00   280918.28     81.0                     0.04
    ## 101    74.00   511081.31     90.0                        -
    ## 102    83.00   371954.66     87.0                        -
    ## 103    93.00   197464.94     78.0                        -
    ## 104    53.00  2759878.73     86.0                     0.06
    ## 105    86.00  3516494.97     88.0                 0.000654
    ## 106    79.00  1043494.11     95.0                        -
    ## 107    39.00  9308720.70     95.0                     0.03
    ## 108    89.00 81303641.95    100.0                  0.00229
    ## 109    79.00 10118357.04     95.0                     0.01
    ## 110    73.00 19049781.56     99.0                     0.03
    ## 111    89.00   775496.78     92.0                     0.07
    ## 112    97.00  2125641.40     99.0                        -
    ##     comblib_sidney2008_rank netmhcpan_ic50 netmhcpan_rank
    ## 1                       0.3              -              -
    ## 2                       7.3              -              -
    ## 3                       4.6              -              -
    ## 4                       0.7              -              -
    ## 5                       2.2              -              -
    ## 6                       1.1              -              -
    ## 7                         -              -              -
    ## 8                        22              -              -
    ## 9                       2.1              -              -
    ## 10                      1.1              -              -
    ## 11                       29              -              -
    ## 12                        1              -              -
    ## 13                       22              -              -
    ## 14                       21              -              -
    ## 15                      6.4              -              -
    ## 16                       75              -              -
    ## 17                       39              -              -
    ## 18                      0.3              -              -
    ## 19                      8.7              -              -
    ## 20                      5.6              -              -
    ## 21                        -              -              -
    ## 22                       16              -              -
    ## 23                       38              -              -
    ## 24                       38              -              -
    ## 25                       53              -              -
    ## 26                      1.4              -              -
    ## 27                       31              -              -
    ## 28                       18              -              -
    ## 29                       67              -              -
    ## 30                       66              -              -
    ## 31                       75              -              -
    ## 32                        -              -              -
    ## 33                       67              -              -
    ## 34                       41              -              -
    ## 35                       22              -              -
    ## 36                        -              -              -
    ## 37                       25              -              -
    ## 38                       32              -              -
    ## 39                        -              -              -
    ## 40                       78              -              -
    ## 41                       50              -              -
    ## 42                       20              -              -
    ## 43                        -              -              -
    ## 44                       46              -              -
    ## 45                       16              -              -
    ## 46                        -              -              -
    ## 47                       35              -              -
    ## 48                       36              -              -
    ## 49                       14              -              -
    ## 50                       29              -              -
    ## 51                       78              -              -
    ## 52                       36              -              -
    ## 53                       41              -              -
    ## 54                       74              -              -
    ## 55                       79              -              -
    ## 56                       43              -              -
    ## 57                        -              -              -
    ## 58                       77              -              -
    ## 59                        -              -              -
    ## 60                       14              -              -
    ## 61                       44              -              -
    ## 62                       80              -              -
    ## 63                       85              -              -
    ## 64                       37              -              -
    ## 65                       24              -              -
    ## 66                        -              -              -
    ## 67                       36              -              -
    ## 68                        -              -              -
    ## 69                       29              -              -
    ## 70                       63              -              -
    ## 71                       86              -              -
    ## 72                        -              -              -
    ## 73                        -              -              -
    ## 74                       58              -              -
    ## 75                       30              -              -
    ## 76                      5.4              -              -
    ## 77                        -              -              -
    ## 78                       87              -              -
    ## 79                       24              -              -
    ## 80                       59              -              -
    ## 81                        -              -              -
    ## 82                       19              -              -
    ## 83                        -              -              -
    ## 84                       61              -              -
    ## 85                        -              -              -
    ## 86                       18              -              -
    ## 87                       34              -              -
    ## 88                        -              -              -
    ## 89                       64              -              -
    ## 90                        -              -              -
    ## 91                       73              -              -
    ## 92                       83              -              -
    ## 93                        -              -              -
    ## 94                       93              -              -
    ## 95                       37              -              -
    ## 96                       70              -              -
    ## 97                        -              -              -
    ## 98                        -              -              -
    ## 99                        -              -              -
    ## 100                      82              -              -
    ## 101                       -              -              -
    ## 102                       -              -              -
    ## 103                       -              -              -
    ## 104                      98              -              -
    ## 105                      27              -              -
    ## 106                       -              -              -
    ## 107                      88              -              -
    ## 108                      53              -              -
    ## 109                      90              -              -
    ## 110                      91              -              -
    ## 111                      96              -              -
    ## 112                       -              -              -
