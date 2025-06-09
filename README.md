# crisprHAL — Better data for better predictions: data curation improves deep learning for sgRNA/Cas9 prediction
Note: this is a static paper-specific repository, and as such, may not contain up-to-date models.

## Links to the crisprHAL repositories, papers, and web tool:
* [Up-to-date crisprHAL prediction tool for use](https://github.com/tbrowne5/crisprHAL)
* [Online crisprHAL prediction tool](https://crisprhal.streamlit.app/)
* [crisprHAL 2.0 paper repository](https://github.com/tbrowne5/Better-data-for-better-predictions-data-curation-improves-deep-learning-for-sgRNA-Cas9-prediction/) — **YOU ARE HERE**
* crisprHAL 2.0 pre-print (Available soon)
* crisprHAL SaCas9 paper repository (Available soon)
* crisprHAL SaCas9 pre-print (Available soon)
* [crisprHAL SpCas9 paper repository](https://github.com/tbrowne5/A-generalizable-Cas9-sgRNA-prediction-model-using-machine-transfer-learning)
* [crisprHAL SpCas9 publication](https://doi.org/10.1038/s41467-023-41143-7)

## ABSTRACT

The Cas9 enzyme along with a single guide RNA molecule is a modular tool for genetic engineering and has shown effectiveness as a species-specific anti-microbial. The ability to accurately predict on-target cleavage is critical as activity varies by target. Using the sgRNA nucleotide sequence and an activity score, predictive models have been developed with the best performance resulting from deep learning architectures. Prior work has emphasized robust and novel architectures to improve predictive performance. Here, we explore the impact of a data-centric approach through optimization of the input target site adjacent nucleotide sequence length and the use of data filtering for read counts in the control conditions  to improve input data utility. Using the existing crisprHAL architecture, we develop crisprHAL Tev, a best-in-class bacterial SpCas9 prediction model with performance that generalizes between related species and across data types. During this process, we also rebuild two prior \textit{E. coli} Cas9 datasets, demonstrating the importance of data quality, and resulting in the production of an improved bacterial eSpCas9 prediction model.

## FILES AND DIRECTORIES:
* data/ contains raw data, model training data, and hold-out test sets
* figures/ contains the R scripts and data for paper figure generation
* models/ contains the model saves
* crisprHAL.py contains the main Python script for model running
* model.py contains the code for running the models
* processing.py contains the code for input data processing

## QUICK START

If you wish to run the model on your own nucleotide sequence follow parts 0 to 3. 

If you wish to validate the model or to update the model with your own data, follow parts 4 to 5.

**Please be advised:** This is a paper-specific repository, for practical use with the most up-to-date models, please visit the [crisprHAL](https://github.com/tbrowne5/crisprHAL) repository.

## Sections of this guide:

Setting up and running the model to predict sgRNA activities:
* 0: Model requirements ```Time: 1-10 minutes```
* 1: Running the model test ```Runtime: 10 seconds```
* 2: Processing nucleotide sequences into model inputs ```Runtime: Varible```
* 3: Predicting with the model ```Runtime: 1-10 seconds```

Additional information and methods: 
* 4: Preparing your own model input files & comparing predictions
* 5: Validating the trained models ```Variable runtime```
* 6: Data availability and processing
* 7: Citations
* 8: How to cite this model

## 0: Requirements

These are in a file called requirements.txt and should be in the working directory.
```
python
tensorflow
```

These can be instantiated within a conda environment: ```Time: 1-10 minutes```

```
conda create --name HAL
conda activate HAL
conda install --file requirements.txt
```

This installation has been tested in Ubuntu 20.04.4 and Mac OSX 10.14.5, but has not been tested on Windows.

## 1: Run model test
```
python crisprHAL.py
```
Test our primary SpCas9/TevSpCas9 prediction model: crisprHAL Tev

Success here is that the model runs without error, showing that it is installed correctly. ```Runtime: ~10 seconds```

## 2: Understand options

```
python crisprHAL.py [options]
```
```
--model, -m   [TevSpCas9, eSpCas9, WT-SpCas9]  Specify the model name (default: TevSpCas9)
```
Model (default=TevSpCas9): specify the model to be used. TevSpCas9 should be used for both TevSpCas9 and SpCas9 predictions. WT-SpCas9 should only be used for crisprHAL WT validation.
```
--input, -i   [Input file path]                Input file for prediction (fasta, csv, or tsv)
```
Input: crisprHAL accepts three types of input files: FASTA (.fasta and .fa), CSV (.csv), and TSV (.tsv). If no input is specified, the model will default to testing on its respective hold-out set.

```
--output, -o  [Output file path]               Output file for prediction results
```
Output: specify the output path and file name of choice. If no output is specified, output file will have the input file name with the prefix: "_predictions.txt"
```
--circular                                     Process fasta as a circular input sequence
```
Circular (default=FALSE): specific to FASTA inputs; specifies that the input sequence should be treated as circular DNA rather than linear.
```
--compare, -c                                  Compare predictions with scores in the input file second column
```
Compare (default=FALSE): specific to CSV/TSV inputs; specifies that the input file contains a second column with scores for comparison. Outputs Spearman and Pearson correlation coefficients between predictions and provided scores, and writes both the predictions and scores to the output file.
```
--train, -t                                    Train the model specified
```
Train (default=FALSE): train the model specified by ```--model/-m [modelName]``` (default=TevSpCas9) using the corresponding training dataset in the data/ directory.
```
--epochs, -e  [Integer epoch value]            Specify number of epochs for training (default: model-specific)
```
Epochs: specify the number of training epochs to be run when training the model. By default each model uses its respective 5CV-identified epochs.
```
--help, -h                                     Show this help message
```
Help: prints available options.

## 3: Predict with the model


## 4: Preparing your own input CSV Files

Input CSV file for prediction only, no "Compare" option:
```
ATGCATATCCCTCTTATTGCCGGTCGCG
GTCTTTATCAGCTAACCAGTCGGTATCC
CGATGGTCAATATCAGCCGTTGGCGCAG
TTTGCCTCATCAACACCTGAAGGCCTCA
CCGCTTTTCCTGCCTGACCTTGGGTGAA
CATAATAGTATTTCCGATAAGGGTCCCC
CGTAGCCGACAATACGGAAACGGTGAGT
CCGAGACGTTGATGCCAATATGGAAATT
CAGGAAACGGCTAACAGAACCGGACCAA
GTGGCAATCGTCGTTTTAACCGGCAAAC
```

Input CSV file for prediction and comparison of the predictions to scores in column 2:
```
CTCGATTGAGGGGCTGGGAATGGGTGAT,8.21062788839
ATCTTTATCGGTCTTAGCAAAGGCTTTG,30.1092205446
CGGGCCAGACTGGCCGAGACGGGTCGTT,11.0586722001
CAGCATCATGCCGTGATCGCCGGGAAAG,0.67958058668
GGGGAAAGGACAAGCAGTGCGGGCAAAA,29.2338752707
GAGAGAATTTTGACCTGATCCGGCTCGC,2.28311187737
CGTGATGCAACTGTGTAATGCGGCTGAC,27.8665335936
GAACATCACCGCCTCACGTCCGGTTTTG,26.3480104838
TCGATTGAGGGGCTGGGAATGGGTGATC,41.2590972746
CCGTGTAAGGGAGATTACACAGGCTAAG,4.25926295656
```

## 5: Validate the training of the model

This will assess whether the training model is working. It will not change the model used for predictions.

Perform 5-fold cross validation with the TevSpCas9 dataset transfer learning from the eSpCas9 base model: ```Runtime: ~20 seconds```
```
python crisprHAL.py train TevSpCas9
```

Perform 5-fold cross validation with the SpCas9 dataset transfer learning from the eSpCas9 base model: ```Runtime: ~20 seconds```
```
python crisprHAL.py train SpCas9
```

Perform 80:20 train-test split with the eSpCas9 dataset (Guo et al. 2018) used as the base model\*: ```Runtime: 1-10 minutes```
```
python crisprHAL.py train eSpCas9
```
\*An 80% training & 20% test split was used for base model generation, and therefore has been included in place of the 5-fold cross validation tests used for the TevSpCas9 and SpCas9 enzyme transfer learning-based models.

## 6: Data availability and processing

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699

NCBI SRA Bioproject: PRJNA939699

Information about data processing can be found under crisprHAL/data/processing.


## 7: Citations

* Guo, J. et al. Improved sgRNA design in bacteria via genome-wide activity profiling. Nucleic acids research **46**, 7052–7069 (2018).
* Abadi, M. et al. TensorFlow: Large-scale machine learning on heterogeneous systems (2015). Software available from tensorflow.org.
* Fernandes, A. D. et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome **2**, 1–13 (2014).
* Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in python. Nat. Methods **17**, 261–272 (2020).

## 8: How to cite crisprHAL

Ham, D.T., Browne, T.S., Banglorewala, P.N. et al. A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets. **Nat Commun** *14*, 5514 (2023). https://doi.org/10.1038/s41467-023-41143-7
