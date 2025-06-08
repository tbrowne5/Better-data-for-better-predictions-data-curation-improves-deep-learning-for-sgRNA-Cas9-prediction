
import numpy as np
from scipy.stats import spearmanr, pearsonr

class processing:
    
    def __init__(self):
        pass

    modelVersionInputLength = {
        "TEVSPCAS9": [37, 3, 14],
        "WT-SPCAS9": [400, 190, 190],
        "ESPCAS9": [300, 140, 140]
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

    def find_targets(self, sequence, inputParameters, circular=True):
        
        if len(sequence) < inputParameters[0]:
            print(f"Error: Input sequence '{sequence}' is shorter than the required length of {inputParameters[0]} bases.")
            return []
        else:
            sequence = sequence.upper().replace("U", "T").replace(" ", "").replace("\n", "")
            if circular: sequence = sequence + sequence[:inputParameters[0] - 1]
            targets = []
            # Find NGG targets
            for i in range(20+inputParameters[1], len(sequence) - inputParameters[2] + 1):
                if sequence[i+1:i+3] == "GG":
                    target = sequence[i-(20+inputParameters[1]):i+inputParameters[2]]
                    if len(target) == inputParameters[0]:
                        targets.append(target)
                    else:
                        print(len(target))
                        break
            return targets

    def process_fasta(self, fastaFile, inputParameters=[37, 3, 14], circular=True):
        
        # Read the FASTA file
        inputSequences = []
        with open(fastaFile, 'r') as file:
            sequence = ""
            for line in file:
                if line.startswith('>'):
                    if len(sequence) > 0:
                        inputSequences += self.find_targets(sequence, inputParameters, circular)
                        # ALSO NEED TO ADD IN THE REVERSE COMPLEMENT OF THE SEQUENCE
                        sequence = ""
                else:
                    sequence += line.strip()
            inputSequences += self.find_targets(sequence, inputParameters, circular)
        
        print(f"Found {len(inputSequences)} target sequences in the FASTA file.")
        return inputSequences, self.onehotencode(inputSequences), None

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

    def read_training_data(self, modelName): return self.process_csv(self, self.modelVersionTrainingData[modelName], predictionColumn=True)

    def read_testing_data(self, modelName): return self.process_csv(self, self.modelVersionTestingData[modelName], predictionColumn=True)

    def read_input(self, modelName, inputFile, compare, circular=True):

        # If no parameters specified, the program will provide hold-out testing data performance
        if inputFile == None:
            return self.process_csv(self.modelVersionTestingData[modelName], predictionColumn=True)
        # Recommended input format for model use (no fiddling with input sequence lengths)
        elif inputFile.endswith('.fasta') or inputFile.endswith('.fa'):
            return self.process_fasta(inputFile, self.modelVersionInputLength[modelName], circular)
        # Method for validating model performance, requires explicit input sequence lengths
        elif inputFile.endswith('.csv') or inputFile.endswith('.tsv'):
            return self.process_csv(inputFile, predictionColumn=compare)
        else:
            raise ValueError("Unsupported file format, please provide an input of the format: .fasta, .fa, .csv, or .tsv")
    
    def compare_predictions(self, predictions, actual):
        
        if len(predictions) != len(actual):
            raise ValueError("Predictions and actual scores must have the same length.")
        
        # Calculate Spearman and Pearson correlations from numpy arrays
        spearman_corr, _ = spearmanr(predictions, np.array(actual))
        pearson_corr, _ = pearsonr(predictions, np.array(actual))

        print(f"Spearman correlation: {spearman_corr:.4f}")
        print(f"Pearson correlation: {pearson_corr:.4f}")
