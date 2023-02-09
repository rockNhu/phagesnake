#!/bin/env python3
#-* coding = UTF-8 *-
# Author = Shixuan Huang (Rock Nhu)
import re,sys
input_file = sys.argv[1]
with open(input_file) as f:
    content = f.read()
res = re.findall('(>(.*?)\n[A-Z|\n|\ ]*)',content)
if res == []:
    exit('The fasta file is empty.')
# simplify the name and save as dict
name_total_name = {}
name_seq = {}
for result in res:
    seq = result[0]
    total_name = result[1]
    name = total_name.split(' ',1)[0]
    name_total_name[name] = total_name
    name_seq[name] = seq.replace(total_name,name)
output_totalname = input_file.rsplit('.',1)[0] + '_totalname.pydict'
with open(output_totalname,'w') as f:
    f.write(str(name_total_name))
output_nameseq = input_file.rsplit('.',1)[0] +'_nameseq.pydict'
with open(output_nameseq,'w') as f:
    f.write(str(name_seq))