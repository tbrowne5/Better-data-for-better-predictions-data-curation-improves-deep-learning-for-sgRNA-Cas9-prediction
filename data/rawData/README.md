
## Generation of scores with target site adjacent nucleotides at increasing minimum control read count values

**Warning:** running the listed R scripts will generate 100 dataset versions with a minimum control read count cutoff from 1–to–100.

All datasets were generated using: R version 4.3.2

### The Data Available Here:
* *C. rodentium* TevSpCas9 sequencing count file with corresponding sgRNA targets
* *C. rodentium* TevSpCas9 target sites with 20 nucleotides upstream and 38 nucleotides (PAM+35) downstream
* *E. coli* eSpCas9 sequencing counts with corresponding sgRNA targets
* *E. coli* eSpCas9 target sites with 540 nucleotides upstream and 558 nucleotides downstream
* *E. coli* WT-SpCas9 sequencing counts with corresponding sgRNA targets
* *E. coli* WT-SpCas9 target sites with 540 nucleotides upstream and 558 nucleotides downstream

### Required R Packages:
* Bioconductor
* ALDEx2
* reshape2
