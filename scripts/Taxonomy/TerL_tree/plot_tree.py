import re
from Bio import Phylo
from matplotlib import pyplot as plt
import sys
tree = Phylo.read(sys.argv[1], "newick")
tree.rooted = True
name_num = len(re.findall(' name=',str(tree)))
fig = plt.figure(figsize=(10, name_num / 3.8))
ax1 = plt.subplot(111)
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['right'].set_visible(False)
Phylo.draw(tree,axes=ax1,do_show=False,ylabel=('',),xlabel=('',),yticks=([],))
plt.savefig(sys.argv[2],dpi=600)
