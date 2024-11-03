import os
import argparse
parser = argparse.ArgumentParser(description='This script is used to ')
parser.add_argument('-i','--input',required=True,help='Path of input dir')
parser.add_argument('-o','--output',default='.',required=False,help='Path to output, default is .')
args = parser.parse_args()
count = 0
out = ''
for file in os.listdir(args.input):
    filename = file.strip().split('_')[0]
    out += filename + ':\n'
    for line in open(os.path.join(args.input,file)):
        if line.startswith('#'):
            continue
        elif not line:
            continue
        count += 1
        out += line
print(f"Total drug resistance gene and viral factor number:{count}")
with open(args.output, 'w') as f:
    f.write(out)