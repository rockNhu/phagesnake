import argparse
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
'''
usage website:
https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer
'''


class MyCustomTranslator(BiopythonTranslator):
    """Custom translator implementing the following theme:

    - Color terminators in green, CDS in blue, all other features in gold.
    - Do not display features that are restriction sites unless they are BamHI
    - Do not display labels for restriction sites
    - For CDS labels just write "CDS here" instead of the name of the gene.

    """

    def to_color(self, anno):
        '''color the annotation'''
        # add the color reference below here, 2:DNA 4:package 5:structure 10:lysin 3:support
        color1 = [
            'hypo', 'Hypo', 'unknow', 'Unknow', 'Uncharactrize', 'family', 'unannotated'
        ]
        color2 = [
            'deaminase', 'kinase', 'polymerase', 'helicase', 'transferase',
            'oxidoreductase', 'pyrophosphatase', 'exonuclease', 'endonuclease',
            'binding', 'ejectosome', 'anti-repressor', 'endodeoxyribonuclease',
            'recombination', 'primase', 'regulator', 'methylase',
            'DNA replication', 'NrdH-redoxin', 'recombinase', 'polymerase',
            'replisome', 'transcription', 'HNH', 'Regulator'
        ]
        color5 = [
            'terminase', 'Terminase', 'packaging', 'protease', 'assembly',
            'joining', 'adaptor', 'DNA maturase'
        ]
        color4 = [
            'scaffold', 'capsid', 'core protein', 'inner', 'collar', 'particle',
            'portal', 'tape', 'virion', 'prohead', 'head', 'tail', 'connector',
            'structural', 'neck', 'Connector'
        ]
        color10 = [
            'lysin', 'holin', 'glycosyl hydrolase', 'lysogentic', 'amidase',
            'lytic', 'N-acetyltransferase', 'peptidase', 'M15A', 'glutaredoxin',
            'Glutaredoxin', 'Lys', 'glycerophosphodiester phosphodiesterase'
        ]
        color3 = [
            'membrane', 'N-myc-interactor', 'Ig-like', 'zinc finger', 'host',
            'integrase', 'Integrase', 'ammonium transporter', 'translocase',
            'transposase', 'Cytoplasm'
        ]
        if any(c in anno for c in color1):
            return 'grey'
        elif any(c in anno for c in color2):
            return 'red'
        elif any(c in anno for c in color10):
            return 'orange'
        elif any(c in anno for c in color4):
            return 'skyblue'
        elif any(c in anno for c in color5):
            return 'blue'
        elif any(c in anno for c in color3):
            return 'green'
        else:
            return 'grey'

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return self.to_color(feature.qualifiers['product'][0])
        elif feature.type == "tRNA":
            return "yellow"
        else:
            return BiopythonTranslator.compute_feature_color(self, feature)

    def compute_feature_label(self, feature):
        unknow_names = ['hypo', 'unannotated']
        if feature.type == "CDS" and any(i in feature.qualifiers['product'][0] for i in unknow_names):
            return None  # do not display hypothetical protein
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if feature.type not in ["gene", 'source']
        ]


class arrow_plot(object):
    '''This class is for doing what.'''

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='This script is used to plot arrow CDS genome')
        parser.add_argument(
            '-i', '--input', required=True,
            help='Path of input gbk'
        )
        parser.add_argument(
            '--gc',
            required=False, default='on',
            help='Switch the gc% plot on or off'
        )
        parser.add_argument(
            '--out_type',
            required=False, default='line',
            help='Switch the output is line or circular'
        )
        parser.add_argument(
            '-o', '--output',
            required=False, default='custom_bopython_translator.pdf',
            help='Path to output, default is custom_bopython_translator.png'
        )
        self.args = parser.parse_args()
        self.size = self.size_controller()
        self.main()

    def size_controller(self):
        '''Using the genome size to set plot width'''
        handle = SeqIO.read(open(self.args.input),'genbank')
        seq_len = len(handle.seq)
        return seq_len * 0.0006
    
    def plot_arrow_n_gc(self):
        '''Plot the local gc content '''
        def gc_plot(record, windows=50):
            '''we use 50bp windows'''
            def gc(s): return 100.0 * \
                len([c for c in s if c in "GC"]) / windows
            xx = np.arange(len(record.sequence) - windows)
            yy = [gc(record.sequence[x: x + windows]) for x in xx]
            return xx, yy

        graphic_record = MyCustomTranslator().translate_record(self.args.input)
        _, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(self.size, 6), sharex=True,
            gridspec_kw={"height_ratios": [4, 1]}
        )
        graphic_record.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)
        xx, yy = gc_plot(graphic_record)
        ax2.fill_between(xx + 25, yy, alpha=0.3)
        ax2.set_ylim(bottom=0)
        ax2.set_ylabel("GC(%)")
        plt.savefig(self.args.output)

    def plot_arrow_in_cir(self):
        '''Plot and export a circular view of the constract'''
        graphic_record = MyCustomTranslator().translate_record(
            self.args.input, record_class=CircularGraphicRecord
        )
        ax2, _ = graphic_record.plot(figure_width=4)
        ax2.figure.tight_layout()
        ax2.figure.savefig(
            "graphic_record_defined_by_hand_circular.png",
            bbox_inches="tight"
        )

    def plot_arrow_in_line(self):
        graphic_record = MyCustomTranslator().translate_record(self.args.input)
        ax, _ = graphic_record.plot(figure_width=self.size)
        ax.figure.tight_layout()
        ax.figure.savefig(self.args.output)

    def main(self):
        if self.args.out_type == 'circular':
            self.plot_arrow_in_cir()
        elif self.args.gc == 'on':
            # gc and cir are unavailable temporarily, because cir seems not good
            self.plot_arrow_n_gc()
        else:
            self.plot_arrow_in_line()


if __name__ == '__main__':
    arrow_plot()
