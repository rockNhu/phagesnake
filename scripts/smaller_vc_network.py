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
    once = []
    for line in open(nwk):
        content = line.strip('\n').split('\t')
        data1 = content[0]
        data2 = content[1]
        if data1 == sample:
            once.append(data2)
        elif data2 == sample:
            once.append(data1)
    twice = set()
    twice_line = []
    for line in open(nwk):
        content = line.strip('\n').split('\t')
        data1 = content[0]
        data2 = content[1]
        if data1 in once or data2 in once:
            twice.add(data1)
            twice.add(data2)
            twice_line.append(line)
    twice_ovv = [
        line for line in open(ovv)
        if any(i in line.split(',')[0] for i in twice)
    ]
    twice_db = [
        line for line in open(db)
        if any(i in line.split('\t')[0] for i in twice)
    ]
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