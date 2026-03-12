# Class 19: Cancer Mutation
Yane Lee PID A17670350

I downloaded my sequence form the class website and moved it tto my
working directory:

From a **BLASTp** search against human refseq I see that it is:

-Official Symbol: TP53 -Official Full Name: tumor protein 53

``` r
library(bio3d)

a <- read.fasta("A17670350_mutant_seq.fa")
read.fasta("A17670350_mutant_seq.fa")
```

                   1        .         .         .         .         .         60 
    wt_healthy     MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    mutant_tumor   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    mutant_tumor   DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYVGSYGERLGFLHSGTAK
                   ******************************************* **** *********** 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    mutant_tumor   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    mutant_tumor   RCSDSDGLAPPRHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
                   *********** ************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    mutant_tumor   YCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
                    *********************************************************** 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    mutant_tumor   PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .  393 
    wt_healthy     GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    mutant_tumor   GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
                   ********************************* 
                 361        .         .         .  393 

    Call:
      read.fasta(file = "A17670350_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 393 position columns (393 non-gap, 0 gap) 

    + attr: id, ali, call

We could score residue conservation

``` r
mutation.sites <- which(conserv(a) < 1)
```

``` r
paste(a$ali[1, mutation.sites],
mutation.sites,
a$ali[2, mutation.sites])
```

    [1] "Q 104 V" "F 109 E" "Q 192 R" "S 241 Y"
