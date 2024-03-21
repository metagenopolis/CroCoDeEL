import sys
import argparse
import pandas as pd
import numpy as np

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--mgs_profiles_path", type=str, help="MGS profiles path.", required=True)
    parser.add_argument("-o", "--output_directory_path", type=str, help="Output directory path.", required=True)
    arguments = parser.parse_args(args)
    return arguments

def main():
    arguments = get_arguments()
    
    mgs_profiles = pd.read_csv(arguments.mgs_profiles_path, sep='\t', header=0, index_col=0)
    number_of_samples = len(mgs_profiles.columns.tolist())
    
    contamination_results = pd.read_csv(arguments.output_directory_path + "/contamination_results.txt", sep='\t')
    
    summary = contamination_results["rate"].describe()

    number_of_contaminations = int(summary["count"])
    minimum = np.round(summary["min"], 4)
    first_quantile = np.round(summary["25%"], 4)
    mean = np.round(summary["mean"], 4)
    third_quantile = np.round(summary["75%"], 4)
    maximum = np.round(summary["max"], 4)
    
    number_of_contaminated_samples = len(contamination_results["target"].unique())
    
    template = "CroCoDeEL - summary report\n\nNumber_of_samples: {number_of_samples}\nNumber_of_contaminations: {number_of_contaminations}\n" \
        "Number_of_contaminated_samples: {number_of_contaminated_samples}\n\nContamination_rate:\nmin: {minimum}\nfirst_quantile: {first_quantile}\n" \
            "mean: {mean}\nthird_quantile: {third_quantile}\nmax: {maximum}"

    infos = {
            'number_of_samples': number_of_samples,
            'number_of_contaminations': number_of_contaminations,
            'number_of_contaminated_samples': number_of_contaminated_samples,
            'minimum': minimum,
            'first_quantile': first_quantile,
            'mean': mean,
            'third_quantile': third_quantile,
            'maximum': maximum
        }

    file = template.format(**infos)
    
    with open('summary_report.ini', 'w') as f:
        f.write(file)
    
    
if __name__ == '__main__':
    main()