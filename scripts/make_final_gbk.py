import pathlib
import re
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
import argparse


class make_final_gbk(object):
    '''This class is for doing what.'''

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='This script is to make the final gbk file using faa, fasta with blastp and eggnog.')
        parser.add_argument('-a', '--faa', required=True,
                            help='Path of input faa')
        parser.add_argument('-f', '--fasta', required=True,
                            help='Path of input fasta')
        parser.add_argument('-b', '--blastp', required=True,
                            help='Path of input blastp_fmt file')
        parser.add_argument('-e', '--eggnog', required=True,
                            help='Path of input eggnog ')
        parser.add_argument('-o', '--output', required=True,
                            help='Path of output gbk ')
        self.args = parser.parse_args()
        self.main()

    def _check_unknow(self, anno):
        unknow_list = ['Uncharacterized','hypo','unknow','family','-','putative']
        return any(un in anno for un in unknow_list)


    def parse_blastp_anno(self, blastp_out):
        name_data = {}
        for line in open(blastp_out):
            content = line.strip('\n').split('\t')
            name = content[0]
            product = content[1]
            if 'Uncharacterized' in product or 'unknown' in product or 'hypothetical' in product:
                continue
            elif name not in name_data:
                name_data[name] = product
        return name_data

    def parse_eggnog_anno(self, eggnog_out):
        name_data = {}
        for line in open(eggnog_out):
            if '#' in line:
                continue
            content = line.strip('\n').split('\t')
            name = content[0]
            product = content[7]
            if (
                'Uncharacterized' in product
                or 'hypothetical' in product
                or 'unknown' in product
                or product == '-'
            ):
                product = 'hypothetical protein'
            gos = content[9]
            ec = content[10]
            kos = content[11]
            cazy = content[-3]
            pfam = content[-1]
            name_data[name] = [product, ec, gos, kos, cazy, pfam]
        return name_data

    def get_faa(self, faa_files):
        '''get orf reigion and fmt the records'''
        name_data = {}
        faa = pathlib.Path(faa_files).read_text()
        res = re.findall('>(.*?)\n([A-Z|\n]*)', faa)
        for result in res:
            header = result[0].split(' # ')
            name = header[0]
            name_data[name] = SeqFeature()
            name_data[name].type = 'CDS'
            start = int(header[1]) - 1
            end = int(header[2])
            strand = int(header[3])
            name_data[name].location = FeatureLocation(start, end, strand)
            desc = header[4].split(';')
            name_data[name].qualifiers['start_type'] = desc[2].split('=')[1]
            name_data[name].qualifiers['rbs_motic'] = desc[3].split('=')[1]
            name_data[name].qualifiers['rbs_spacer'] = desc[4].split('=')[1]
            name_data[name].qualifiers['product'] = 'hypothetical protein'
            name_data[name].qualifiers['translation'] = result[1].replace(
                '\n', '')
        return name_data

    def new_record(self, sample_name, fasta_path, faa_data, eggnog_data, blastp_data):
        '''make the final gbk records'''
        new_record = []
        prot_id = 1
        for record in SeqIO.parse(fasta_path, 'fasta'):
            seq_id = record.id
            features = record.features
            for name, data in faa_data.items():
                # {name:[product, ec, gos, kos, cazy, pfam]}
                eggnog = eggnog_data.get(name)
                blastp = blastp_data.get(name)  # product
                contig = name.rsplit('_', 1)[0]  # the contig name
                if seq_id == contig:
                    feature = data
                    feature.qualifiers['locus_tag'] = '{}_{:04d}'.format(
                        sample_name, prot_id)
                    prot_id += 1
                    if eggnog is not None:
                        self._extract_eggnog_feature(eggnog, feature)
                    if blastp is not None:  # finally using the blastp annotation
                        feature.qualifiers['product'] = blastp
                    feature.qualifiers['transl_table'] = '11'
                    features.append(feature)
            record.annotations['molecule_type'] = "DNA"
            new_record.append(SeqRecord(record.seq, seq_id, record.name,
                record.description, record.dbxrefs, features, record.annotations))
        return new_record

    def _extract_eggnog_feature(self, eggnog, feature):
        feature.qualifiers['product'] = eggnog[0]  # description
        if eggnog[1] != '-':
            feature.qualifiers['EC'] = eggnog[1]
        if eggnog[2] != '-':
            feature.qualifiers['GO'] = eggnog[2]
        if eggnog[3] != '-':
            feature.qualifiers['KEGG_ko'] = eggnog[3]
        if eggnog[4] != '-':
            feature.qualifiers['CAZy'] = eggnog[4]
        if eggnog[5] != '-':
            feature.qualifiers['Pfam'] = eggnog[5]

    def main(self):
        faa = self.args.faa
        fasta_path = self.args.fasta
        sample_name = faa.rsplit('.', 1)[0].rsplit('/', 1)[1]
        faa_data = self.get_faa(faa)
        eggnog_anno = self.parse_eggnog_anno(self.args.eggnog)
        blastp_anno = self.parse_blastp_anno(self.args.blastp)
        new_records = self.new_record(
            sample_name, fasta_path, faa_data, eggnog_anno, blastp_anno)
        SeqIO.write(new_records, self.args.output, 'genbank')


if __name__ == '__main__':
    make_final_gbk()
