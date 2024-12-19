#!/bin/env python3
#-* coding = UTF-8 *-
# @Author = Shixuan Huang
import os
import pathlib
import argparse

parser = argparse.ArgumentParser(
    description='This script is used to get neibour nucleotide and fmt the names')
parser.add_argument('-i', '--input', required=True,
                    help='the output dir: blastn_out')
parser.add_argument('-o', '--output', required=True, help='combined fasta')
args = parser.parse_args()

def fmt_name(name):
    if ',' in name:
        name = name.split(',')[0]
    name = name.replace(' ', '_').replace('UNVERIFIED:_', '')
    return name

gate = args.input
out = ''
for file in os.listdir(gate):
    if file.endswith('.fasta'):
        head_num = 1
        for line in open(os.path.join(gate, file)):
            if head_num:
                head_num -= 1
                name = line.strip('>\n')
                newname = fmt_name(name)
                out += '>' + newname + '\n'
                continue
            out += line
with open(args.output, 'w') as f:
    f.write(out)

