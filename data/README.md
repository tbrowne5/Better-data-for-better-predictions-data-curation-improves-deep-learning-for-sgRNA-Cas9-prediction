
## Model Data and Raw Data

This directory contains the crisprHAL Tev, crisprHAL eSp, and crisprHAL WT model training and hold-out test sets, with model-specific input lengths.

### Files in the (TEV/ESP/WT) subdirectory:
* pTox TevSpCas9 independent test set with requisite model input length
* pTox WT-SpCas9 independent test set with requisite model input length
* KatG TevSpCas9 independent test set with requisite model input length

### Files in the RawData subdirectory:
* Experimental sequence count files used for score generation
* sgRNA/Cas9 target sites with adjacent nucleotide sequences
* R scripts for re-generating all control condition read count curation datasets — Warning: the R scripts will generate 100 datasets each

### Files in the SupplementaryData subdirectory:
* S1: Entire datasets for A) TevSpCas9, B) eSpCas9, C) WT-SpCas9
* S2: Long target sites for A) TevSpCas9, B) eSpCas9 and WT-SpCas9
* S3: Sequencing read count tables for the 3 primary datasets
* S4: ALDEx2 ```aldex.effect``` outputs for the 3 datasets
* S5: Independent datasets with crisprHAL Tev input lengths
* S6: Sequencing read count tables for the 3 independent datasets

