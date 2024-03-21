from datetime import datetime
import sys
import argparse

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--mgs_profiles_path", type=str, help="MGS profiles path.", required=True)
    parser.add_argument("-o", "--output_directory_path", type=str, help="Output directory path.", required=True)
    arguments = parser.parse_args(args)
    return arguments

def main():
    arguments = get_arguments()
    
    template = "CroCoDeEL - v. 1.0.0\n\n[Execution]\ndate: {date}\n\n[Options]\nMGS_profiles: {mgs_profiles}\n\n" \
    "[Ouput]\nresults: {contamination_results}\nresults_plot: {contamination_results_plot}\n"

    date = datetime.now().strftime("%Y-%m-%d")

    infos = {
            'date': date,
            'mgs_profiles': arguments.mgs_profiles_path,
            'contamination_results': arguments.output_directory_path + "/contamination_results.txt",
            'contamination_results_plot': arguments.output_directory_path + "/contamination_results.pdf",
        }

    file = template.format(**infos)
    
    with open('parameters.ini', 'w') as f:
        f.write(file)
    
    
if __name__ == '__main__':
    main()