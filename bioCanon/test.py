import pandas as pd
import os
from bioCanon.utils.vcf_helper import sanitize_vcf_file

class pandas_vcf:
    required_fields = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    def find_header_line(self, file):
        line_num = 0
        with open(file) as search:
            for line in search:
                if line[0:6] == '#CHROM':
                    return [line_num,line.rstrip()]
                line_num += 1
        return [-1,'']

    def find_missing_fields(self):
        header = self.get_header()
        missing_fields = []
        for field in self.required_fields:
            if field not in header:
                missing_fields.append(field)
        return missing_fields

    def set_header(self,header):
        self.header = header

    def get_header(self):
        return self.header

    def set_sample_list(self,samples):
        self.sample_list = samples

    def get_sample_list(self):
        return self.sample_list

    def check_file(self, file):
        if not os.path.isfile(file):
            return False
        if os.path.getsize(file) == 0:
            return False
        return True

    def check_basecall_format(self):
        sample_id = self.get_header()[0]
        return self.df[sample_id].dtype

    def set_positions(self,positions):
        self.positions = positions

    def get_positions(self):
        return self.positions

    def get_df(self):
        return self.df

    def process(self):
        positions = {}
        base_lookup = {}
        base_callformat = self.check_basecall_format()

        convert = False
        if base_callformat == 'int64':
            convert = True

        self.df = self.df.astype(str)
        self.df['POS'] = self.df['POS'].astype('int64')
        for index,row in self.df.iterrows():
            chr = row['#CHROM']
            pos = row['POS']
            ref_base = row['REF']
            alt_bases = row['ALT'].split(',')

            if not chr in positions:
                positions[chr] = []
                base_lookup[chr] = {}
            positions[chr].append(pos)
            base_lookup[chr][pos] = []


            #convert numeric bases into nucleotides
            bases = [ref_base]
            for base in alt_bases:
                if base not in ['A','T','C','G']:
                    base = 'N'
                bases.append(base)

            base_lookup[chr][pos] = bases

        self.set_positions(positions)



    def __init__(self,file):

        if not self.check_file(file):
            return None

        (line_num,header_line) = self.find_header_line(file)
        self.set_header(header_line.split("\t"))
        missing_fields = self.find_missing_fields()
        if len(missing_fields) > 0:
            return None
        self.set_sample_list(list(set(self.get_header()) - set(self.required_fields)))
        if line_num == 0:
            skip_lines = 0
        else:
            skip_lines = line_num


        self.df = pd.read_csv(file,skip_blank_lines=True, header=line_num,skiprows=skip_lines,sep="\t")
        self.df.columns = self.get_header()
        self.process()




pandas_vcf("/Users/jrobertson/PycharmProjects/bioCanon/bioCanon/data/snp.sites.gaps.txt").check_basecall_format()