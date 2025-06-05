import sys
from model import model
from processing import process

training = False
modelName = "TEVSPCAS9"
modelNames = ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9", "WTSPCAS9", "WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT", "ESPCAS9", "ESP"]
epochs = 50
inputFile = None
outputFile = None
compare = False


# PERHAPS MOVE THIS TO PROCESSING AND HAVE A VERY SIMPLE CRISPRHAL.PY FILE!!!

def parse_args(args):
    global training, modelName, modelNames, epochs, inputFile, outputFile, compare

    for i in range(len(args)):
        if args[i] == "--train" or args[i] == "-t": training = True
        elif args[i] == "--input" or args[i] == "-i": inputFile = args[i + 1]
        elif args[i] == "--output" or args[i] == "-o": outputFile = args[i + 1]
        elif args[i] == "--compare" or args[i] == "-c": compare = args[i + 1]
        elif args[i] == "--epochs" or args[i] == "-e": epochs = int(args[i + 1])
        elif args[i] == "--model" or args[i] == "-m" or args[i] == "--enzyme":
            if args[i + 1].upper() in modelNames:
                if args[i + 1].upper() in ["TEVSPCAS9","TEVCAS9", "TEV", "SPCAS9"]: modelName = "TEVSPCAS9"
                if args[i + 1].upper() in ["WT-SPCAS9", "WILD-TYPE-SPCAS9","WILD-TYPE","WILDTYPE", "WT"]: modelName = "WT-SPCAS9"
                if args[i + 1].upper() in ["ESPCAS9", "ESP"]: modelName = "ESPCAS9"
            else:
                print(f"Error: Model '{args[i + 1]}' is not recognized. Available models: TevSpCas9, eSpCas9, and WT-SpCas9")
                sys.exit(1)
        elif args[i] == "--help" or args[i] == "-h":
            print("Usage: python crisprHAL.py [options]")
            print("Options:")
            print("  --model, -m         Specify the model name (default: TevSpCas9)")
            print("  --input, -i         Input file for training or prediction")
            print("  --output, -o        Output file for results")
            print("  --compare, -c       Input file contains second column of scores for comparison")
            print("  --train, -t         Train the model")
            print("  --epochs, -e        Number of training epochs (default: 50)")
            print("  --help, -h          Show this help message")
            sys.exit(0)
    
    if training == True and inputFile is not None:
        print("Error: Please specify either training mode or an input file for prediction generation.")
        sys.exit(1)

# Read in arguments
# If no input specified, use default hold-out test set
# If no output specified, but input specified, strip file type annd add _predictions.csv

def run_model():
    global training, modelName, inputFile, outputFile, compareFile, epochs
    
    model.load_model(modelName + ".keras")

    # Model name provides input sequence length for processing
    # If inputFile default of "None" is passed, the hold-out test set will be used instead
    # The compare flag indicates that the input file contains a second column of scores to be used for comparison
    predictionData = process.read_input(modelName, inputFile, compare)
    
    if training:
        print("Training model")
    else:
        predictionData = model.predict(predictionData)

if __name__ == "__main__":
    parse_args(sys.argv[1:])
    run_model()