from plotnine import ggplot, aes, geom_point, geom_abline, coord_cartesian, theme, element_text, labs, scale_x_log10, scale_y_log10, theme_minimal, theme_void, scale_color_manual, guides
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
import patchworklib as pw
import sys
import argparse

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--input", type=str, help="MGS profiles path.", required=True)
    #parser.add_argument("-c", "--contamination_results", type=str, help="File of detected contaminations.", required=True)
    parser.add_argument("-o", "--output_directory_path", type=str, help="Output directory path.", required=True)
    arguments = parser.parse_args(args)
    return arguments

def main():
    arguments = get_arguments()
    
    mgs_profiles = pd.read_csv(arguments.input, sep='\t', header=0, index_col=0)
    mgs_profiles = mgs_profiles.div(mgs_profiles.sum(axis=0), axis=1)

    contamination_results = pd.read_csv(arguments.output_directory_path + "/contamination_results.txt", sep='\t')
    contamination_results = contamination_results.sort_values(by=['probability', 'rate'], ascending=[False, False])

    nb_plots = contamination_results.shape[0]
    nb_plots_to_add = (16 - (nb_plots % 16)) % 16

    plots = []
    size = 10

    if contamination_results.shape[0] > 0:
        for i in range(contamination_results.shape[0]):
            sample1 = str(contamination_results.iloc[i, :]["source"])
            sample2 = str(contamination_results.iloc[i, :]["target"])
            probability = float(contamination_results.iloc[i, :]["probability"])
            intercept = -np.log10(float(contamination_results.iloc[i, :]["rate"]))
            contamination_rate = round(float(contamination_results.iloc[i, :]["rate"]) * 100, 2)
            msp_list = contamination_results.iloc[i, :]["msp_list"]
            msp_list = msp_list.split(',')

            data = mgs_profiles[[sample1, sample2]]
            data.columns = ['source', 'target']

            spearman = data.corr(method='spearman').iloc[0, 1]

            min_value = data.loc[(data['target'] != 0) & (data['source'] != 0)].min().min()
            max_value = data.loc[(data['target'] != 0) & (data['source'] != 0)].max().max()

            data['category'] = ['contamination_line' if index in msp_list else 'other_point' for index in data.index]
            
            g = (ggplot(data) +
                geom_point(aes(x='target', y='source', color='category'), size=1, shape='o', fill='none') +
                geom_abline(intercept=0, slope=1, color="grey", size=0.1) +
                geom_abline(aes(intercept=intercept, slope=1), linetype="dashed", color="red", size=0.1, alpha=0.5) +
                coord_cartesian(xlim=(np.log10(min_value/10), np.log10(max_value*10)),
                                ylim=(np.log10(min_value/10), np.log10(max_value*10))) +
                scale_x_log10() +
                scale_y_log10() +
                scale_color_manual(values={'contamination_line': 'orange', 'other_point': 'black'}) +
                theme_minimal() +
                theme(axis_title_x=element_text(size=size),
                    axis_title_y=element_text(size=size),
                    axis_text_x=element_text(size=size),
                    axis_text_y=element_text(size=size),
                    title=element_text(size=size)) +
                labs(x=sample2, 
                    y=sample1,
                    title="prob = {}, rate = {}%, rho = {}".format(probability, contamination_rate, round(spearman, 2))) +
                guides(color=False))
            
            pw_g = pw.load_ggplot(g, figsize=(3,3))

            plots.append(pw_g)
            
    for i in range(nb_plots_to_add):
        g = (ggplot() + theme_void())
        pw_g = pw.load_ggplot(g, figsize=(3,3))
        plots.append(pw_g)
    
    p = PdfPages("contamination_results.pdf") 
    for page in range(len(plots)//16):
        rows = []
        for k in range(page*16, page*16+16, 4):
            row = pw.stack(plots[k:k+4], operator="|", margin=0.2)
            rows.append(row)
        all_page = pw.stack(rows, operator='/', margin=0.2)
        all_page.savefig(p, format='pdf')
    p.close()
    
if __name__ == '__main__':
    main()