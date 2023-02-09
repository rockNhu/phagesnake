import argparse
from dna_features_viewer import BiopythonTranslator
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
        #add the color reference below here, 2:DNA 4:package 5:structure 10:lysin 3:support
        color1 = ['hypo','Hypo','unknow','Unknow','Uncharactrize','family']
        color2 = ['deaminase','kinase','polymerase','helicase','transferase','oxidoreductase',
            'pyrophosphatase','exonuclease','endonuclease','binding','ejectosome','anti-repressor',
            'endodeoxyribonuclease','recombination','primase','regulator','methylase',
            'DNA replication','NrdH-redoxin','recombinase','polymerase','replisome',
            'transcription','HNH','Regulator']
        color5 = ['terminase','Terminase','packaging','protease','assembly','joining','adaptor','DNA maturase']
        color4 = ['scaffold','capsid','core protein','inner','collar','particle','portal','tape','virion','prohead',
            'head','tail','connector','structural','neck','Connector']
        color10 = ['lysin','holin','glycosyl hydrolase','lysogentic','amidase','lytic','N-acetyltransferase',
            'peptidase','M15A','glutaredoxin','Glutaredoxin','Lys','glycerophosphodiester phosphodiesterase']
        color3 = ['membrane','N-myc-interactor','Ig-like','zinc finger','host','integrase',
            'Integrase','ammonium transporter','translocase','transposase','Cytoplasm']
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
        if feature.type == "CDS" and 'hypo' in feature.qualifiers['product'][0]:
            return None # do not display hypothetical protein
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        """Do not display promoters. Just because."""
        return [
            feature for feature in features
            if feature.type not in ["gene", 'source']
        ]

parser = argparse.ArgumentParser(description='This script is used to plot arrow CDS genome')
parser.add_argument('-i','--input',required=True,help='Path of input gbk')
parser.add_argument('-o','--output',default='custom_bopython_translator.pdf',required=False,help='Path to output, default is custom_bopython_translator.png')
args = parser.parse_args()
graphic_record = MyCustomTranslator().translate_record(args.input)
ax, _ = graphic_record.plot(figure_width=20)
ax.figure.tight_layout()
ax.figure.savefig(args.output)