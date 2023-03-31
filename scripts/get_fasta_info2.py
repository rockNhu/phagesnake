#!/bin/env python3
# -* coding = UTF-8 *-
# @Author = Shixuan Huang (Rock Nhu)
import os
import argparse


class get_GC(object):
    '''This class is for counting GC % & length.'''

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='This script is used to count GC % & length & scaffolds_num')
        parser.add_argument('-i', '--input', required=True,
                            help='Path of input dir')
        parser.add_argument('-orf', '--orf_nums', default='NA',
                            required=False, help='Path of input orf.list')
        parser.add_argument('-o', '--output', default='.',
                            required=False, help='Path to output, default is .')
        args = parser.parse_args()
        self.input_gate = args.input
        if args.orf_nums != 'NA':
            self.orfs = self.get_orf(args.orf_nums)
            self.head = 'Name\tGC%\tLength\tORFs\n'
        else:
            self.orfs = None
            self.head = 'Name\tGC%\tLength\n'
        self.na_ph = self.get_input()
        self.main(args.output)

    def get_input(self):
        na_ph = []
        for file in os.listdir(self.input_gate):
            name = os.path.basename(file).rsplit('.', 1)[0]
            path = os.path.join(self.input_gate, file)
            na_ph.append([name, path])
        return na_ph

    def get_orf(self, orf_list):
        orfs = {}
        for line in open(orf_list):
            content = line.strip().split(':')
            name = content[0]
            num = content[1]
            orfs[name] = num
        return orfs

    def count(self, seq):
        total = len(seq)
        if not total:
            return False, False
        n = sum(i in ['g', 'c', 'G', 'C'] for i in seq)
        gc = n / total * 100
        return gc, total

    def main(self, output):
        output_line = []
        for bsname, path in self.na_ph:
            the_seq = ''.join(line.strip()
                              for line in open(path) if '>' not in line)
            gc, length = self.count(the_seq)
            if not gc:
                if bsname:
                    print(f"Error: divied by 0:{bsname}")
                continue
            if self.orfs is None:
                output_line.append(f'{bsname}\t{gc}\t{length}\n')
            else:
                orf = self.orfs[bsname]
                output_line.append(f'{bsname}\t{gc}\t{length}\t{orf}\n')
        with open(f'{output}/seq_info.tsv', 'w') as f:
            f.write(self.head)
            f.writelines(output_line)


if __name__ == '__main__':
    get_GC()
