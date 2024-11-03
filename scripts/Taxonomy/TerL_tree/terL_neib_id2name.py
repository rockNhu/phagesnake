#!/bin/env python3
# -* coding = UTF-8 *-
# @Author = Shixuan Huang
import argparse
import pathlib
import re


class id2name(object):
    '''This class is for doing what.'''

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='This script is used to catch genome name,\
            to specify the phage with their scientific names')
        parser.add_argument('-i', '--input', required=True,
                            help='Path of input faa')
        parser.add_argument('-tn', '--totalname', required=True,
                            help='Path of input genome totalname.pydict')
        parser.add_argument('-o', '--output', required=True,
                            help='Path of output')
        self.args = parser.parse_args()
        self.ref = self.get_ref()
        self.main()

    def get_ref(self):
        '''get name dict'''
        content = pathlib.Path(self.args.totalname).read_text()
        return eval(content)

    def main(self):
        content = pathlib.Path(self.args.input).read_text()
        out = ''
        result = re.findall('>(.*?)\n([A-Z|\n]*)', content)
        for res in result:
            contig_name = res[0]
            # make sure the seq is endswith '\n'
            seq = res[1] if res[1].endswith('\n') else res[1] + '\n'
            name = contig_name.split('_', 1)[0]
            # translate to name and remove id and space
            new_name = self.ref[name].split(' ',1)[1].replace(' ','_')
            # shorten name, remove right with ','
            new_name = new_name.rsplit(',',1)[0]
            out += f'>{new_name}\n{seq}'
        with open(self.args.output, 'w') as o:
            o.write(out)


if __name__ == '__main__':
    id2name()
