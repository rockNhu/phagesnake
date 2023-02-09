#!/bin/env python3
#-* coding = UTF-8 *-
# @Author = Shixuan Huang
import os
import argparse
class get_seqs_from_dict(object):
    '''This class is for doing what.'''
    def __init__(self):
        parser = argparse.ArgumentParser(description='This script is used to catch fasta and seperate them to a dir or all of them to a file')
        parser.add_argument('-i','--input',required=True,help='Path of input list')
        parser.add_argument('-ns','--nameseqs',required=True,help='Path of input nameseqs.pydict')
        parser.add_argument('-tn','--totalname',required=True,help='Path of input totalname.pydict')
        parser.add_argument('-o','--output',required=True,help='Path of output')
        parser.add_argument('--seperate',help='output to a dir or a file')
        args = parser.parse_args()
        self.main(args)


    def extract_fasta_sep(self,nameseqs,totalname,to_ex,output):
        if not os.path.exists(output):
            os.makedirs(output)
        for i in to_ex: # i is id
            seq = nameseqs[i]
            name = totalname[i].replace(' ','_').replace('/','_')
            out = seq.replace(i,name)
            with open(os.path.join(output, f'{name}.fasta'), 'w') as f:
                f.write(out)

    def extract_fasta_one(self,nameseqs,totalname,to_ex,output):
        with open(output, 'w') as f:
            for i in to_ex: # i is id
                seq = nameseqs[i]
                name = totalname[i].replace(' ','_').replace('/','_')
                out = seq.replace(i,name)
                f.write(out)

    def exclude_empty(self,to_ex,args):
        if not to_ex:
            print(f'Empty input: {args.input}')
            if args.seperate:
                if not os.path.exists(args.output):
                    os.makedirs(args.output)
            else:
                with open(args.output,'w') as f:
                    f.write('')

    def main(self,args):
        # load_data
        to_ex = [line.strip() for line in open(args.input)]
        self.exclude_empty(to_ex,args)
        with open(args.nameseqs, 'r') as f:
            nameseqs = eval(f.read()) # {id:>id\nseqs\n}
        with open(args.totalname, 'r') as f:
            totalname = eval(f.read())# {id:totalname}
        # main
        if args.seperate:
            self.extract_fasta_sep(nameseqs,totalname,to_ex,args.output) # this output is dir
        else:
            self.extract_fasta_one(nameseqs,totalname,to_ex,args.output) # this output is file

if __name__ == '__main__':
    get_seqs_from_dict()
