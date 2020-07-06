from version import __version__
import logging
from argparse import ArgumentParser
from bioCanon.version import __version__


LOG_FORMAT = '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'

def init_console_logger(lvl):
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description='BIO_Hansel Scheme developement')

    parser.add_argument('--in_vcf', type=str, required=True, help='VCF File of SNPs')
    parser.add_argument('--in_nwk', type=str, required=False, help='Newick Tree of strains')
    parser.add_argument('--in_meta', type=str, required=False, help='Tab delimited file of metadata', default=None)
    parser.add_argument('--in_groups', type=str, required=False, help='Tab delimited file of (sample_id|genotype)', default=None)
    parser.add_argument('--reference', type=str, required=True, help='Reference fasta sequence from VCF')
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--prefix', type=str, required=False, help='Prefix for output files',default='bioCanon')
    parser.add_argument('--root_name', type=str, required=False, help='Name of sample to root tree',default='')
    parser.add_argument('--root_method', type=str, required=False, help='Method to root tree (midpoint,outgroup)', default=None)
    parser.add_argument('--max_missing_position', type=float, required=False, help='Maximum percentage of samples which can be missing a position for it to be valid',
                        default=0.25)
    parser.add_argument('--max_missing_sample', type=float, required=False, help='Maximum percentage of positions which a sample can be missing for it to be included',
                        default=0.25)
    parser.add_argument('--min_polymorphism', type=float, required=False, help='Minimum number of samples with different states for a position to be valid',
                        default=0.25)
    parser.add_argument('--fconst', type=str, required=False, help='Number of constant sites suitable for IQ-TREE -fconst',
                        default=None)
    parser.add_argument('--num_threads', type=int, required=False,
                        help='Number of threads to use',
                        default=1)
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()



def main():
    cmd_args = parse_args()



logging = init_console_logger(3)
vcf_file = '/Users/jrobertson/Google Drive/__SH_Subtyping_Manuscript/__2020-06-HeidelbergManuscript/testing/core.full.vcf'
filtered_vcf = '/Users/jrobertson/Google Drive/__SH_Subtyping_Manuscript/__2020-06-HeidelbergManuscript/testing/test.vcf'
tree_file = '/Users/jrobertson/Google Drive/__SH_Subtyping_Manuscript/__2020-06-HeidelbergManuscript/testing/raxml.final_tree.tre'

tree = parse_tree(tree_file, logging,set_root=True,resolve_polytomy=True,method='outgroup',outgroup='ERR230402_snippy')
sanitize_vcf_file(vcf_file,out_file=filtered_vcf,snp_log_file='/Users/jrobertson/Desktop/test_snps.log',sample_log_file='/Users/jrobertson/Desktop/test_samples.log',logging=logging,min_count=1)
#vcf_to_fasta('/Users/jrobertson/Google Drive/__SH_Subtyping_Manuscript/__2020-06-HeidelbergManuscript/core.vcf','/Users/jrobertson/Desktop/test.fasta',position_filter=[],sample_filter=[])
#clades = check_cansnp_support(filtered_vcf, tree, logging,homoplasy_report_file='/Users/jrobertson/Google Drive/__SH_Subtyping_Manuscript/__2020-06-HeidelbergManuscript/testing/homoplasy.txt')
