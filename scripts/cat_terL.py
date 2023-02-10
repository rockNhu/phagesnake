#!/bin/env python3
# -* coding = UTF-8 *-
# @Author = Shixuan Huang
import re
import argparse
import pathlib


class get_terL(object):
    '''This class is to cat TerL from self TerL and neibour TerL,
    output is final TerLs faa.'''

    def __init__(self):
        self.get_input_var()
        self.main()
        self.output()

    def get_input_var(self):
        parser = argparse.ArgumentParser(description='Format name & Cat')
        parser.add_argument('--self', required=True,
                            help='Path of input self terl faa')
        parser.add_argument('--neib', required=True,
                            help='Path of input neibour terl faa')
        parser.add_argument('-o', '--output', required=True,
                            help='Path of output faa')
        parser.add_argument('-s', '--sample', required=True, help='name')
        self.args = parser.parse_args()
        # parse the input variable
        self.phage_terLs = {}  # TODO: merge the split terL seq
        self.new_content = []  # TODO: the neibour TerL contnet to output
        self.seq = '' # TODO: the self TerL contnet to output

    def main(self):
        # fmt output contig name
        content = pathlib.Path(self.args.neib).read_text()
        res = re.findall('>(.*?)\n([A-Z|\n]*)', content)
        for result in res:
            name = result[0]
            seq = result[1].replace('\n', '')
            if len(seq) < 150:  # remove fragment
                continue
            phage_name = name.split(' ')[1] if ' ' in name else name
            # remove unknown phage
            if '_sp.' in phage_name or 'uncultured' in phage_name:
                continue
            # concat multiple terL to a single
            if phage_name not in self.phage_terLs:
                self.phage_terLs[phage_name] = seq
            else:
                self.phage_terLs[phage_name] += seq
        # neibour output
        self.new_content = [f'>{phage_name}\n{seq}\n' \
                for phage_name, seq in self.phage_terLs.items()]
        # get self output
        self.seq = pathlib.Path(self.args.self).read_text()

    def output(self):
        # write the final output
        with open(self.args.output, 'w') as f:
            # the neibours TerL contnet
            f.writelines(self.new_content)
            # the self TerL content
            if self.seq != '':
                f.write(f'>{self.args.sample}\n{self.seq}\n')
            else:
                f.write('')


if __name__ == '__main__':
    get_terL()
