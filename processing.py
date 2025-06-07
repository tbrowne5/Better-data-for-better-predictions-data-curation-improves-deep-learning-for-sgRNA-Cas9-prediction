
import numpy as np

modelVersionInputLength = {
    "TEVSPCAS9": [37, 3, 14],
    "ESPCAS9": [400, 190, 190],
    "WT-SPCAS9": [300, 140, 140]
    }

modelVersionPath = {
    "TEVSPCAS9": "models/TevSpCas9.h5",
    "ESPCAS9": "models/eSpCas9.h5",
    "WT-SPCAS9": "models/WT-SpCas9.h5"
    }

modelVersionTrainingData = {
    "TEVSPCAS9": "data/TevSpCas9_training_data.csv",
    "ESPCAS9": "data/eSpCas9_training_data.csv",
    "WT-SPCAS9": "data/WT-SpCas9_training_data.csv"
    }

modelVersionTestingData = {
    "TEVSPCAS9": "data/TevSpCas9_testing_data.csv",
    "ESPCAS9": "data/eSpCas9_testing_data.csv",
    "WT-SPCAS9": "data/WT-SpCas9_testing_data.csv"
    }

modelVersionDefaultEpochs = {
    "TEVSPCAS9": 50,
    "ESPCAS9": 50,
    "WT-SPCAS9": 50
    }

def onehotencode(self, nucleotideEncoding):
    
    nt_to_ohe_dict = {"A": [1, 0, 0, 0], "C": [0, 1, 0, 0], "G": [0, 0, 1, 0], "T": [0, 0, 0, 1], "M": [0, 0, 0, 0], "N": [0, 0, 0, 0]}
    onehotencoding = []
    for sequence in nucleotideEncoding: onehotencoding.append(np.array([nt_to_ohe_dict[nt.upper()] for nt in sequence]))
    return np.array(onehotencoding)

def process_fasta(self, fastaFile, inputParameters=[37, 3, 14]):
    return

def process_csv(self, csvFile, predictionColumn=False):
    # Check if comma or tab separated via file name
    if csvFile.endswith('.csv'): delimiter = ','
    elif csvFile.endswith('.tsv'): delimiter = '\t'
    else: raise ValueError("Unsupported file format. Please provide a .csv or .tsv file.")
    
    # Read the CSV file
    inputSequences = []
    scores = []
    
    with open(csvFile, 'r') as file:
        for line in file:
            if line.strip():
                targetAndScore = line.strip().split(delimiter)
                inputSequences.append(targetAndScore[0])
                if predictionColumn and len(targetAndScore) > 1: scores.append(float(targetAndScore[1]))
                else: scores.append(targetAndScore[0])

    oneHotEncoded = self.onehotencode(inputSequences)

    if predictionColumn: return inputSequences, oneHotEncoded, np.array(scores)
    else: return inputSequences, oneHotEncoded, None

def read_training_data(self, modelName): return process_csv(self, modelVersionTrainingData[modelName], predictionColumn=True)

def read_testing_data(self, modelName): return process_csv(self, modelVersionTestingData[modelName], predictionColumn=True)

def read_input(self, modelName, inputFile, compare):

    # If no parameters specified, the program will provide hold-out testing data performance
    if inputFile == None:
        return self.process_csv(modelVersionTestingData[modelName], predictionColumn=True)
    # Recommended input format for model use (no fiddling with input sequence lengths)
    elif inputFile.endswith('.fasta') or inputFile.endswith('.fa'):
        return self.process_fasta(inputFile)
    # Method for validating model performance, requires explicit input sequence lengths
    elif inputFile.endswith('.csv') or inputFile.endswith('.tsv'):
        return self.process_csv(inputFile, predictionColumn=compare)
    else:
        raise ValueError("Unsupported file format, please provide an input of the format: .fasta, .fa, .csv, or .tsv")
    
