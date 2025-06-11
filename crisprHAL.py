import sys
import os

os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "3" # Specifies which GPU to use
os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async" # may need: export TF_GPU_ALLOCATOR=cuda_malloc_async

from models import models, modelVersionDefaultEpochs
from processing import processing

training = False
modelName = "TEVSPCAS9"
modelNames = ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9", "WTSPCAS9", "WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT", "ESPCAS9", "ESP"]
epochs = None
inputFile = None
outputFile = None
compare = False
circularInput = False
summary = False
replicate = False
replicationSets = ["pTox_TevSpCas9.csv", "pTox_WTSpCas9.csv", "KatG_TevSpCas9.csv"]

# PERHAPS MOVE THIS TO PROCESSING AND HAVE A VERY SIMPLE CRISPRHAL.PY FILE!!!

def parse_args(args):
    global training, modelName, modelNames, epochs, inputFile, outputFile, compare, circularInput, summary, replicate

    for i in range(len(args)):
        if args[i] == "--train" or args[i] == "-t": training = True
        elif args[i] == "--input" or args[i] == "-i": inputFile = args[i + 1]
        elif args[i] == "--output" or args[i] == "-o": outputFile = args[i + 1]
        elif args[i] == "--circular": circularInput = True
        elif args[i] == "--compare" or args[i] == "-c": compare = True
        elif args[i] == "--epochs" or args[i] == "-e": epochs = int(args[i + 1])
        elif args[i] == "--summary" or args[i] == "-s": summary = True
        elif args[i] == "--replicate" or args[i] == "-r": replicate = True
        elif args[i] == "--model" or args[i] == "-m" or args[i] == "--enzyme":
            if args[i + 1].upper() in modelNames:
                if args[i + 1].upper() in ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9"]: modelName = "TEVSPCAS9"
                if args[i + 1].upper() in ["WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT"]: modelName = "WT-SPCAS9"
                if args[i + 1].upper() in ["ESPCAS9", "ESP"]: modelName = "ESPCAS9"
            else:
                print(f"Error: Model '{args[i + 1]}' is not recognized. Available models: TevSpCas9, eSpCas9, and WT-SpCas9.\nFor up-to-date and other models including SaCas9 please use: github.com/tbrowne5/crisprHAL")
                sys.exit(1)
        elif args[i] == "--help" or args[i] == "-h":
            print("Usage: python crisprHAL.py [options]")
            print("Options:")
            print("  --model, -m   [TevSpCas9, eSpCas9, WT-SpCas9]  Specify the model name (default: TevSpCas9)")
            print("  --input, -i   [Input file path]                Input file for prediction (fasta, csv, or tsv)")
            print("  --output, -o  [Output file path]               Output file for prediction results")
            print("  --circular                                     Process fasta as a circular input sequence")
            print("  --compare, -c                                  Compare predictions with scores in the input file second column")
            print("  --train, -t                                    Train the model specified")
            print("  --epochs, -e  [Integer epoch value]            Specify number of epochs for training (default: model-specific)")
            print("  --help, -h                                     Show this help message")
            sys.exit(0)
    
    if training == True and inputFile is not None:
        print("Error: Please specify either training mode or an input file for prediction generation.")
        sys.exit(1)

# Read in arguments
# If no input specified, use default hold-out test set
# If no output specified, but input specified, strip file type annd add _predictions.csv

def run_model():
    global training, modelName, inputFile, outputFile, compare, epochs, circularInput, summary, replicate, replicationSets, summary

    process = processing()
    model = models(modelName, summary)
    
    if training:
        print("Training model")
        trainingData = process.read_training_data(modelName)
        testingData = process.read_testing_data(modelName)
        if epochs is None: epochs = modelVersionDefaultEpochs[modelName]
        for i in range(0,epochs):
            model.train(trainingData[1], trainingData[2], epochs=1, batch_size=1024, verbose=1)
            process.compare_predictions(model.predict(testingData[1]), testingData[2], message=f"Epoch {i+1} on hold-out test set for {modelName}:")
    else:
        # Model name provides input sequence length for processing
        # If inputFile default of "None" is passed, the hold-out test set will be used instead
        # The compare flag indicates that the input file contains a second column of scores to be used for comparison
        if inputFile is None: compare = True
        inputSequences, encodedInputSequences, inputScores = process.read_input(modelName, inputFile, compare, circularInput)
        print(len(inputSequences))
        model.load_model(modelName)
        predictions= model.predict(encodedInputSequences)
        if compare:
            # Compare predictions with the second column of scores in the input file
            process.compare_predictions(predictions, inputScores)
        process.write_predictions(inputSequences, predictions, outputFile, inputFile, inputScores)
    
    if replicate:
        # Replicate the model and save it
        model.load_model(modelName)
        compare = True
        for replicationSet in replicationSets:
            inputSequences, encodedInputSequences, inputScores = process.read_input(modelName, "data/"+modelName+"_TEST_SETS/"+replicationSet, compare)
            process.compare_predictions(model.predict(encodedInputSequences), inputScores, message=f"Replicated {modelName} model performance on {replicationSet}:")

if __name__ == "__main__":
    print("\nWelcome to crisprHAL 2.0 — Better data for better predictions: data curation improves deep learning for sgRNA/Cas9 prediction\n\nPlease be aware that this is a static repository specific to the paper, and as such may not contain up-to-to date models.\n\nThe production models can be found at github.com/tbrowne5/crisprHAL or online at crisprHAL.streamlit.app — thank you!")
    parse_args(sys.argv[1:])
    run_model()
