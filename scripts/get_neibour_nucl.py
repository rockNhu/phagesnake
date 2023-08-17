import argparse
import pathlib
parser = argparse.ArgumentParser(
    description='This script is used to get neibour nucleotide infomation')
parser.add_argument('-i', '--input', required=True,
                    help='the output of blastn, blastn.tsv')
parser.add_argument('-o_list', '--output_list',
                    required=True, help='blastn.list')
parser.add_argument('-o_fmt', '--output_fmt_tsv',
                    required=True, help='blastn_fmt.tsv')
parser.add_argument('-nucl_id_db', '--nucleotide_id_database', required=True,
                    help='the name:total_name reference')
parser.add_argument('-idt', '--identities', default=75, type=float,
                    required=False, help='the idt filter value of blastn')
parser.add_argument('-cov', '--coverages', default=75, type=float,
                    required=False, help='the cov filter value of blastn')
args = parser.parse_args()

# get real name
total_genome_name = eval(pathlib.Path(args.nucleotide_id_database).read_text())

# extract
to_extract = {}
for line in open(args.input):
    content = line.strip('\n').split('\t')  # input type is tsv
    acc = content[1]  # the second column is accession number
    idt = float(content[2])  # the third column is fraction identity
    cov = float(content[3])  # the fourth column is fraction coverage
    if acc not in to_extract:
        to_extract[acc] = [[idt], [cov]]
    else:
        to_extract[acc][0].append(idt)  # the idtentity list
        to_extract[acc][1].append(cov)  # the coverage list
# remove idt and cov away from one genus
fmt_out = 'Accession\tidentity%\tcoverage%\n'
final_acc = []
for acc, values in to_extract.items():
    final_idt = sum(values[0]) / len(values[0])  # the idt is avg
    final_cov = sum(values[1]) * 100  # sum up cov and convert it to perc
    if final_idt >= args.identities and final_cov >= args.coverages:
        final_acc.append(acc)
    fmt_out += f'{total_genome_name[acc]}\t{final_idt}\t{final_cov}\n'
# output to name list
with open(args.output_list, 'w') as f:
    f.writelines(i + '\n' for i in final_acc)
# output to fmt blastn
with open(args.output_fmt_tsv, 'w') as f:
    f.write(fmt_out)
