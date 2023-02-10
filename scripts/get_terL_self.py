#!/bin/env python3
# -* coding = UTF-8 *-
# @Author = Shixuan Huang
import argparse
import pathlib
import re


class get_TerL(object):
    '''This class is to get TerL from blastp_out,
    get alignment output to find the TerL,
    output is a list of target protein ids.'''

    def __init__(self):
        self.get_input_var()
        self.main()
        self.output()

    def get_input_var(self):
        # parse the input variable
        parser = argparse.ArgumentParser(
            description='This script is used to get TerL.list from blastp_out')
        parser.add_argument('-i', '--input', required=True,
                            help='Path of input blastp_out file')
        parser.add_argument('-o', '--output', required=True,
                            help='Path of output list')
        parser.add_argument('-f', '--faa', required=True,
                            help='Path of input faa, self faa')
        parser.add_argument('-of', '--out_faa', required=True,
                            help='Path of output faa, self terL faa')
        parser.add_argument('-s', '--sample', required=True,
                            help='name')
        self.args = parser.parse_args()
        # set the several description types of "Terminase large subunit"
        self.TerL_desc = [
            'DNA maturase B', 'DNA packaging',
            'terminase large', 'Terminase large',
            'terminase, large subunit', 'TerL',
            'Large Terminase', 'Terminase'
        ]
        self.to_extract = set()  # TODO: make output acc list
        self.self_prot_tag = set()  # TODO: get self terL tag
        self.seq = ''  # TODO: get self terL seqs

    def main(self):
        for line in open(self.args.input):
            content = line.strip('\n').split('\t')  # the input is tsv
            tag = content[0]  # first column is the tag of self protein
            # second column is the description of subject
            description = content[1]
            acc = content[-1]  # final column is the acc of subject
            for ter in self.TerL_desc:
                if ter in description:  # check the annotation is terminase
                    self.to_extract.add(acc)
                    self.self_prot_tag.add(tag)
        # get self faa
        content = pathlib.Path(self.args.faa).read_text()
        # merge all terminase large to one, 
        # some phage have multi terL, because they were "broken"
        try:
            for spt in self.self_prot_tag:
                re_find_seq = re.findall(f'>{spt} .*?\n([A-Z|\n]*)', content)[0]
                if len(re_find_seq) > 150:  # terminase large always > 150 bp
                    self.seq += re_find_seq.replace('\n', '')
        except:
            self.seq = ''
            print(f'Error: Terminase Not Found in {self.args.sample}')

    def output(self):
        # output neibour ids to list
        with open(self.args.output, 'w') as f:
            if self.to_extract != '':
                # the subject protein
                f.writelines(i + '\n' for i in self.to_extract)
            else:  # the target protein is None
                f.write('')
        # output self TerL to faa
        with open(self.args.out_faa, 'w') as f:
            f.write(self.seq)


if __name__ == '__main__':
    get_TerL()
