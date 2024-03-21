import pandas as pd
import sys
import argparse

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--mgs_profiles_path", type=str, help="MGS profiles path.", required=True)
    arguments = parser.parse_args(args)
    return arguments

def main():
    arguments = get_arguments()
    mgs_profiles = pd.read_csv(arguments.mgs_profiles_path, sep='\t', header=0, index_col=0)
    samples = mgs_profiles.columns.tolist()

    with open('samples.csv', 'w') as fichier:
        fichier.write("sample\n")
        for sample in samples:
            fichier.write(str(sample) + '\n')
        
if __name__ == '__main__':
    main()