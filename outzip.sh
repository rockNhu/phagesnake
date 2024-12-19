# catch important output using zip
for i in `ls output`;do
    if [ -d output/${i} ]
    zip -r output.zip \
    output/2.seq_info* \
    output/${i}/1.annotations \
    output/${i}/3.nucleotide_alignment/fmt_blastn.tsv \
    output/${i}/3.nucleotide_alignment/viridic_out/04_VIRIDIC_out \
    output/${i}/4.TerL_phylogenetic_tree
done
# or catch important output using 7z
# for i in `ls output`;do
#     if [ -d output/${i} ]
#     7z a output.7z \
#     output/2.seq_info* \
#     output/${i}/1.annotations \
#     output/${i}/3.nucleotide_alignment/fmt_blastn.tsv \
#     output/${i}/3.nucleotide_alignment/viridic_out/04_VIRIDIC_out \
#     output/${i}/4.TerL_phylogenetic_tree
# done