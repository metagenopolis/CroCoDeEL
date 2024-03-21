import pandas as pd
import sys
import argparse
import numpy as np

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-o", "--output_directory_path", type=str, help="Output directory path.", required=True)    
    parser.add_argument("-m", "--longitudinal_metadata", type=str, help="Metadata of longitudinal samples.", required=True)
    arguments = parser.parse_args(args)
    return arguments

def main():
    arguments = get_arguments()
    
    contamination_results = pd.read_csv(arguments.output_directory_path + "/contamination_results.txt", sep='\t')
    contamination_results = contamination_results.sort_values(by=['probability', 'rate'], ascending=[False, False])

    longitudinal_samples = pd.read_csv(arguments.longitudinal_metadata, sep='\t')
    longitudinal_samples.columns = ['samples', 'individual']

    temp = pd.merge(contamination_results, longitudinal_samples, how='left', left_on="source", right_on="samples", validate='1:m')
    temp = temp.drop(["samples"], axis=1)
    temp = pd.merge(temp, longitudinal_samples, how='left', left_on="target", right_on="samples", validate='1:m')
    temp = temp.drop(["samples"], axis=1)
    temp = temp.rename(columns={'individual_x': 'individual_source', 'individual_y': 'individual_target'})

    temp["longitudinal"] = ['yes' if pd.notna(x) and pd.notna(y) and x == y else 'no' if pd.notna(x) and pd.notna(y) else 'missing_information' for x, y in zip(temp['individual_source'], temp['individual_target'])]

    temp_sorted = temp.sort_values(by=['longitudinal'], ascending=False)
    temp_sorted = temp_sorted.drop_duplicates(subset=['source', 'target'], keep='first')
    temp_sorted = temp_sorted.sort_values(by=['probability', 'rate'], ascending=[False, False])
    
    new_contamination_results = temp_sorted.drop(["individual_source", "individual_target"], axis=1)
    columns = [col for col in new_contamination_results.columns if col != 'msp_list'] + ['msp_list']
    new_contamination_results = new_contamination_results.reindex(columns=columns)
    new_contamination_results.to_csv("contamination_results.txt", sep='\t', index=False)
    
if __name__ == '__main__':
    main()