import os
import re
import pathlib
import argparse
parser = argparse.ArgumentParser(description='This script is used to split fasta')
parser.add_argument('-i','--input',required=True,help='Path of input fna_dir')
parser.add_argument('-o','--output',required=True,help='Path of output fmt_fna_dir')
args = parser.parse_args()

if not os.path.exists(args.output):
    os.makedirs(args.output)

for file in os.listdir(args.input):
    try:
        filename = file.rsplit('.',1)[0]
        content = pathlib.Path(os.path.join(args.input,file)).read_text()
        result = re.findall('>(.*?)\n([A-Z|a-z|\n]*)',content)
        for idx,res in enumerate(result):
            name = res[0].replace(' ','_')
            seq = res[1].upper().replace('-','').replace(' ','')
            genome_length = len(seq.replace('\n',''))
            if genome_length < 5000:
                print(f'The {filename}_{idx} genome length:{genome_length} is lower than 5000, skipped.')
                continue
            if not seq.endswith('\n'):
                seq += '\n'
            with open(f'{args.output}/{filename}_{idx}.fasta','w') as o:
                o.write(f'>{name}\n{seq}')
    except Exception as e:
        print(f"Error: Skipped input file: {file},\n file is not valid FASTA format!\n{e}\n")