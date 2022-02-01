import logging
from bioCanon.constants import LOG_FORMAT
from Bio import SeqIO

def init_console_logger(lvl):
    """

    Parameters
    ----------
    lvl [int] : Integer of level of logging desired 0,1,2,3

    Returns
    -------

    logging object

    """
    logging_levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    report_lvl = logging_levels[lvl]

    logging.basicConfig(format=LOG_FORMAT, level=report_lvl)
    return logging


def read_fasta_dict(fasta_file):
    """

    :param fasta_file: [str] Path to fasta file to read
    :return: [dict] of sequences indexed by sequence id
    """
    seqs = dict()
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[str(record.id)] = str(record.seq).upper()
    handle.close()
    return seqs

def write_scheme(genotype_scheme,outfile):
    out_string = ''
    for chrom in genotype_scheme:
        for pos in genotype_scheme[chrom]:

            out_string += ">{}:{}-{}\n{}\n".format(chrom,
                                                  pos,
                                                  genotype_scheme[chrom][pos]['genotype'],
                                                  genotype_scheme[chrom][pos]['in_kmer'])
            for out_kmer in genotype_scheme[chrom][pos]['out_kmers']:
                out_string += ">negative{}:{}-{}\n{}\n".format(chrom, pos, genotype_scheme[chrom][pos]['genotype'],out_kmer)

    fh = open(outfile, 'w')
    fh.write(out_string)
    fh.close()

