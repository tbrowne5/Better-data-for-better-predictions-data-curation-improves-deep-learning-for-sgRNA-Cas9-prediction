# crisprHAL — [Better data for better predictions: data curation improves deep learning for sgRNA/Cas9 prediction](https://www.biorxiv.org/content/10.1101/2025.06.24.661356v1.full)
Note: this is a static paper-specific repository, and as such, may not contain up-to-date models. For the crisprHAL prediction tool, please visit the first link below or visit the website at the second link below.

<p align="center">
  <img src="https://github.com/tbrowne5/crisprHAL/blob/main/data/crisprHAL_Logo.png" width="600">
</p>

## Links to the related repositories, papers, and web tool:
* [Up-to-date crisprHAL prediction tool repository for use](https://github.com/tbrowne5/crisprHAL)
* [Online crisprHAL prediction tool](https://crisprhal.streamlit.app/) ([Repository](https://github.com/tbrowne5/crisprHAL_streamlit))
* [crisprHAL 2.0 SpCas9 paper repository](https://github.com/tbrowne5/Better-data-for-better-predictions-data-curation-improves-deep-learning-for-sgRNA-Cas9-prediction/) — **YOU ARE HERE**
* [crisprHAL 2.0 SpCas9 pre-print](https://www.biorxiv.org/content/10.1101/2025.06.24.661356v1.full)
* [crisprHAL SaCas9 paper repository](https://github.com/tbrowne5/Adenine-methylated-PAM-sequences-inhibit-SaCas9-activity)
* [crisprHAL SaCas9 pre-print](https://www.biorxiv.org/content/10.1101/2025.08.13.670096v1)
* [crisprHAL SpCas9 paper repository](https://github.com/tbrowne5/A-generalizable-Cas9-sgRNA-prediction-model-using-machine-transfer-learning)
* [crisprHAL SpCas9 publication](https://doi.org/10.1038/s41467-023-41143-7)

## ABSTRACT

The Cas9 enzyme along with a single guide RNA molecule is a modular tool for genetic engineering and has shown effectiveness as a species-specific anti-microbial. The ability to accurately predict on-target cleavage is critical as activity varies by target. Using the sgRNA nucleotide sequence and an activity score, predictive models have been developed with the best performance resulting from deep learning architectures. Prior work has emphasized robust and novel architectures to improve predictive performance. Here, we explore the impact of a data-centric approach through optimization of the input target site adjacent nucleotide sequence length and the use of data filtering for read counts in the control conditions  to improve input data utility. Using the existing crisprHAL architecture, we develop crisprHAL Tev, a best-in-class bacterial SpCas9 prediction model with performance that generalizes between related species and across data types. During this process, we also rebuild two prior *E. coli* Cas9 datasets, demonstrating the importance of data quality, and resulting in the production of an improved bacterial eSpCas9 prediction model.

## FILES AND DIRECTORIES:
* data/ contains model training & testing data, raw data, and supplementary data
* figures/ contains the R scripts and data for paper figure generation
* models/ contains the model saves
* crisprHAL.py contains the main Python script for model running
* models.py contains the code for running the models
* processing.py contains the code for input data processing

## QUICK START

If you wish to run the model on your own nucleotide sequence follow parts 0 to 3. 

If you wish to validate the model predictions or training, follow parts 4 to 5.

If you wish to simply obtain predictions, you can do so easily through the [crisprHAL website](https://crisprHAL.streamlit.app).

**Please be advised:** This is a paper-specific repository, for practical use with the most up-to-date models, please visit the [crisprHAL](https://github.com/tbrowne5/crisprHAL) repository.

## SECTIONS OF THIS GUIDE:

Setting up and running the model to predict sgRNA/Cas9 activities:
* 0: Model requirements ```Time: 1-10 minutes```
* 1: Running the model test ```Runtime: 10 seconds```
* 2: Understanding available model options
* 3: Predicting with the model ```Runtime: 1-10 seconds```

Additional information and methods: 
* 4: Preparing your own model input files & comparing predictions
* 5: Validating the trained models ```Variable runtime```
* 6: Data availability and processing
* 7: Citations
* 8: How to cite this model

## 0: Requirements

Python and Tensorflow versions used:
```
python==3.12
tensorflow==2.19
```
Create a virtual environment to install Tensorflow and SciPy:
```
python3.12 -m venv HAL
source HAL/bin/activate
pip install tensorflow[and-CUDA]==2.19
pip install scipy==1.15
```

This installation has been tested in Ubuntu 20.04.4 and Mac OSX 10.14.5, but has not been tested on Windows.


## 1: Run the model test
```
python crisprHAL.py
```
Test our primary SpCas9/TevSpCas9 prediction model: crisprHAL Tev

Success here is that the model runs without error, showing that it is installed correctly. ```Runtime: ~10 seconds```


## 2: Understand options

```
python crisprHAL.py [options]
```
Options:
```
--model, -m   [TevSpCas9, eSpCas9, WT-SpCas9]
```
Model (default=TevSpCas9): specify the model to be used. TevSpCas9 should be used for both TevSpCas9 and SpCas9 predictions. WT-SpCas9 should only be used for crisprHAL WT validation.
```
--replicate, -r
```
Replicate: load the independent test sets from ```data/[modelName]_TEST_SETS/``` and replicate the model performance shown in the paper.
```
--input, -i   [Input file path]
```
Input: crisprHAL accepts three types of input files: FASTA (.fasta and .fa), CSV (.csv), and TSV (.tsv). If no input is specified, the model will default to testing on its respective hold-out set.

```
--output, -o  [Output file path]
```
Output: specify the output path and file name of choice. If no output is specified, output file will have the input file name with the prefix: "_predictions.txt"
```
--circular
```
Circular (default=FALSE): specific to FASTA inputs; specifies that the input sequence should be treated as circular DNA rather than linear.
```
--compare, -c
```
Compare (default=FALSE): specific to CSV/TSV inputs; specifies that the input file contains a second column with scores for comparison. Outputs Spearman and Pearson correlation coefficients between predictions and provided scores, and writes both the predictions and scores to the output file.
```
--train, -t
```
Train (default=FALSE): train the model specified by ```--model/-m [modelName]``` (default=TevSpCas9) using the corresponding training dataset in the data/ directory.
```
--epochs, -e  [Integer epoch value]
```
Epochs: specify the number of training epochs to be run when training the model. By default each model uses its respective 5CV-identified epochs. Without the ```--train/-t``` flag, ```--epochs/-e``` will do nothing.
```
--help, -h
```
Help: prints available options.


## 3: Predict with the model

To run a test of the model's predictions, please run the command:
```
python crisprHAL.py -m TevSpCas9 -i sample.fa
```
This will: 1) Identify all sgRNA targets in the sample.fa file, 2) Predict Cas9 activity for each target, and 3) Write the targets and predicted activities to a file called sample_predictions.txt

The full list of options for FASTA input-based model predictions are:
```
python crisprHAL.py --model [TevSpCas9/eSpCas9/WT-SpCas9] --input [input file path] --output [optional output file path] --circular
```


## 4: Preparing your own input CSV Files

Preparing your own input CSV (or TSV) file requires more work, as the input must be the exact length required by the specific model. Different models have different target site adjacent sequence inclusion requirements as follows:

* SpCas9/TevSpCas9: ```37 nucleotides = 3 upstream, 20 sgRNA, 14 downstream (NGG PAM plus 11)```
* eSpCas9: ```406 nucleotides = 193 upstream, 20 sgRNA, 193 downstream (NGG PAM plus 190)```
* WT-SpCas9: ```378 nucleotides = 189 upstream, 20 sgRNA, 169 downstream (NGG PAM plus 166)```

Input CSV file for prediction only, SpCas9/TevSpCas9 model input, no "Compare" option:
```
GCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTA
GATACCAGGATCTTGCCATCCTATGGAACTGCCTCGG
GTAACCATGCATCATCAGGAGTACGGATAAAATGCTT
CTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGG
AATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGG
GTGACGACTGAATCCGGTGAGAATGGCAAAAGCTTAT
AAAACGCCATTAACCTGATGTTCTGGGGAATATAAGG
CTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTT
```

Predict using the command: ```python crisprHAL.py -i [fileName.csv]``` with optional specification of model ```-m [modelName]``` and output path ```-o [outputPath]```

Input CSV file for prediction and prediction vs score comparison in column 2, SpCas9/TevSpCas9 model input:
```
GCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTA,1.17419584308754
GATACCAGGATCTTGCCATCCTATGGAACTGCCTCGG,-0.482944213762817
GTAACCATGCATCATCAGGAGTACGGATAAAATGCTT,0.421119225849831
CTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGG,0.654659434852531
AATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGG,-0.995596435390042
GTGACGACTGAATCCGGTGAGAATGGCAAAAGCTTAT,-2.07685671663022
AAAACGCCATTAACCTGATGTTCTGGGGAATATAAGG,-0.294359200889779
CTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTT,0.92532668013273
```

Predict using the command: ```python crisprHAL.py -i [fileName.csv] -c``` with optional specification of model ```-m [modelName]``` and output path ```-o [outputPath]```


## 5: Validate the training of the model

This will assess whether the training model is working. It will not change the model used for predictions.

Train the model for the model-specific default number of epochs. If no model is specified, crisprHAL Tev (SpCas9/TevSpCas9) will be trained.
```
python crisprHAL.py -t -m [modelName]
```

Traing the model for a specific number of epochs:
```
python crisprHAL.py -t -m [modelName] -e [integerValueForTrainingEpochs]
```

Example: train crisprHAL eSp (eSpCas9) for 50 epochs:
```
python crisprHAL.py -m eSpCas9 -t -e 50
```


## 6: Data availability and processing

Raw sequence reads from which the TevSpCas9 and SpCas9 datasets are derived are available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA939699 and https://www.ncbi.nlm.nih.gov/bioproject/PRJNA450978

NCBI SRA BioProject: PRJNA939699

NCBI SRA BioProject: PRJNA450978

Information about data processing can be found under crisprHAL/data/rawData


## 7: Citations

* Guo, J. et al. Improved sgRNA design in bacteria via genome-wide activity profiling. Nucleic acids research **46**, 7052–7069 (2018).
* Abadi, M. et al. TensorFlow: Large-scale machine learning on heterogeneous systems (2015). Software available from tensorflow.org.
* Fernandes, A. D. et al. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome **2**, 1–13 (2014).
* Virtanen, P. et al. SciPy 1.0: fundamental algorithms for scientific computing in python. Nat. Methods **17**, 261–272 (2020).

## 8: How to cite crisprHAL

Browne, T.S., et al. Better data for better predictions: data curation improves deep learning for sgRNA/Cas9 prediction. **BioRxiv** (Pre-Print). https://doi.org/10.1101/2025.06.24.661356

Ham, D.T., Browne, T.S., Banglorewala, P.N. et al. A generalizable Cas9/sgRNA prediction model using machine transfer learning with small high-quality datasets. **Nat Commun** *14*, 5514 (2023). https://doi.org/10.1038/s41467-023-41143-7
