import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import argparse

class clustermap(object):
    def __init__(self):
        parser = argparse.ArgumentParser(description='This script is used to plot pyANI output file ANIb_percentage_identity.tab')
        parser.add_argument('-i','--input',required=True,help='Path of input tab file')
        parser.add_argument('-o','--output',required=True,help='Path to output svg file')
        parser.add_argument('--show_number',default=False,required=False,help='show number in final plot')
        args = parser.parse_args()
        self.input_file = args.input
        self.output_name = args.output
        self.anno = args.show_number
        self.df = self.get_file()
        self.main()
    
    def fmt_name(self, n):
        if ',' in n:
            n = n.split(',')[0]
        if 'phage_' in n:
            n = n.split('phage_')[1]
        if 'virus_' in n:
            n = n.split('virus_')[1]
        if '_genome_assembly' in n:
            n = n.split('_genome_assembly')[0]
        return n

    def get_file(self):
        df = pd.read_csv(self.input_file,sep='\t',index_col=0)
        return df.rename(index=lambda x:self.fmt_name(x),columns=lambda x:self.fmt_name(x))

    def main(self):
        plt.rc('font',family='Times New Roman')
        cmap = mcolors.LinearSegmentedColormap.from_list(
            'custom_map',
            [(0, '#a0a0a0'), (0.6999, '#a0a0a0'), (0.7, 'blue'),(0.9, "blue"), (0.95, 'white'), (1, 'red')]
        )
        if self.anno:
            afont = {
                'family':'Arial',
                'size':10,
                'weight':'normal',
                'color':'black'
                }
            sns.clustermap(self.df,method ='ward',annot=True,fmt='.2',annot_kws=afont,metric='euclidean',cmap=cmap,vmin=0,vmax=1)
        else:
            sns.clustermap(self.df,method ='ward',metric='euclidean',cmap=cmap,vmin=0,vmax=1)
        plt.savefig(self.output_name)

if __name__ == '__main__':
    clustermap()
