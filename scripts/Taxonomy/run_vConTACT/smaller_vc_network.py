#!/bin/env python3
# -* coding = UTF-8 *-
# @Author = Shixuan Huang
import argparse


def parse_args():
    '''Parse all arguments.'''
    parser = argparse.ArgumentParser(
        description='This script is used to easily visualize output of vConTACT2.')
    parser.add_argument('-nwk', '--network', required=True,
                        help='Path of input `c1.ntw`, the main network file.')
    parser.add_argument('-ovv', '--overview', required=True,
                        help='Path of input `genome_by_genome_overview.csv`, the VC and sub_VC infomation.')
    parser.add_argument('-db', '--database', required=True,
                        help='Path of database file `data.tsv`, it conclude total taxo of all phages.')
    parser.add_argument('-s', '--sample', required=True,
                        help='The phage name, the output node will cauculate around it.')
    parser.add_argument('-onwk', '--output_small_ntw', required=True,
                        help='Path of output file.')
    parser.add_argument('-oovv', '--output_small_ovv', required=True,
                        help='Path of output file.')
    parser.add_argument('-odb', '--output_small_db', required=True,
                        help='Path of output file.')
    return parser.parse_args()


def check_singleton(ovv, sample):
    '''Check if the phage was singleton'''
    for line in open(ovv):
        return sample in line and 'Singleton' in line

def filt_twice(nwk, ovv, db, sample):
    '''Get twice neibour of object.'''
    class nwk_reader():
        '''.ntw file has no header'''
        def __init__(self, file):
            self.reader = open(file)
            self.handle = iter(self.reader)

        def __iter__(self):
            for line in self.handle:
                content = line.strip('\n').split(' ')
                data1 = content[0]
                data2 = content[1]
                yield data1, data2, line

        def close(self):
            self.reader.close()


    once = set()
    for data1, data2, _ in nwk_reader(nwk):
        if data1 == sample:
            once.add(data2)
        elif data2 == sample:
            once.add(data1)
    twice = set() # twice neibour is contain once
    for data1, data2, _ in nwk_reader(nwk):
        if data1 in once or data2 in once:
            twice.add(data1)
            twice.add(data2)
    twice_line = [
        line for data1, data2, line in nwk_reader(nwk)
        if data1 in twice and data2 in twice
    ]

    def get_data(file, sp):
        twice_out = []
        head_num = 1
        for line in open(file):
            if head_num:
                head_num -= 1
                twice_out.append(line)
                continue
            if any(i in line.split(sp)[0] for i in twice):
                twice_out.append(line)
        return twice_out
    
    twice_ovv = get_data(ovv, ',')
    twice_db = get_data(db, '\t')

    return twice_line, twice_ovv, twice_db

def main():
    args = parse_args()
    # The singleton is not observed in cytoscape.
    if check_singleton(args.overview, args.sample):
        print('Sample is singleton.')
        nwk_filt = ''
        csv_filt = ''
        metas_filt = ''
    else:
        nwk_filt, csv_filt, metas_filt = filt_twice(args.network, args.overview, args.database, args.sample)
    with open(args.output_small_ntw,'w') as f:
        f.writelines(nwk_filt)
    with open(args.output_small_ovv,'w') as f:
        f.writelines(csv_filt)
    with open(args.output_small_db,'w') as f:
        f.writelines(metas_filt)

if __name__ == '__main__':
    main()