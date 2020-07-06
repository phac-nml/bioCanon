import pandas as pd
import os




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

    def transform_basecalls(self):
        self.df.ix['x', 'C']

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
        print(self.df)


print(pandas_vcf("/Users/jrobertson/PycharmProjects/bioCanon/bioCanon/data/snp.sites.gaps.txt").check_basecall_format())