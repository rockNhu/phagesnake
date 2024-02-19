import argparse
import pathlib


def is_unknown(product):
    unknow_words = [
        'Hypothetical', 'hypothetical','Unknow','unknow', 'unannotated', 'Unannotated'
        'Uncharactrized','uncharactrized','gp','Gp'
    ]
    return any(un in product for un in unknow_words)

parser = argparse.ArgumentParser(
    description='This script is used to get neibour protein infomation')
parser.add_argument('-i', '--input', required=True,
                    help='the output of blastn, blastn.tsv')
parser.add_argument('-o_fmt', '--output_fmt_tsv',
                    required=True, help='blastp_fmt.tsv')
parser.add_argument('-nucl_id_db', '--nucleotide_id_database', required=True,
                    help='the name:total_name reference')
parser.add_argument('-prot_id_db', '--protein_id_database', required=True,
                    help='the protein name:total_name reference')
parser.add_argument('-idt', '--identities', default=70, type=float,
                    required=False, help='the idt filter value of blastn')
parser.add_argument('-cov', '--coverages', default=70, type=float,
                    required=False, help='the cov filter value of blastn')
args = parser.parse_args()

# load the id:totalname
total_genome_name = eval(pathlib.Path(args.nucleotide_id_database).read_text())
total_protein_name = eval(pathlib.Path(args.protein_id_database).read_text())

out_line = 'Name\tDescription\tCoverage%\tIdentity%\tTarget\tTarget_Total_Name\n'
for line in open(args.input):
    contents = line.strip('\n').split('\t')
    name = contents[0]
    target = contents[1]
    product = total_protein_name[target].split(' ',1)[1]
    idt = float(contents[2])
    qcov = float(contents[3])
    target_name = total_genome_name.get(target.split('_')[0], '')
    if not target_name:
        print(f'Error target name: {target}')
        continue
    # the filter
    if not is_unknown(product) and qcov >= args.identities and idt >= args.coverages:
        out_line += f'{name}\t{product}\t{qcov}\t{idt}\t{target}\t{target_name}\n'

with open(args.output_fmt_tsv, 'w') as f:
    f.write(out_line)