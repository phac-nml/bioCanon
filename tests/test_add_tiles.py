from bioCanon import __main__
from Bio import SeqIO
import os


def test_add_tiles():
    conflict_pos = {359: [366], 515: [507], 1895: [1893]}
    scheme = {9: {}, 8: {1: [{'position': 359, 'ref_base': 'T', 'alt_base': 'C', 'ref_count': 17,
                              'alt_count': 5, 'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
                              'qc_warnings': []},
                             {'position': 515, 'ref_base': 'T', 'alt_base': 'G', 'ref_count': 17,
                              'alt_count': 5, 'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
                              'qc_warnings': []},
                             {'position': 9797, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17,
                              'alt_count': 5, 'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
                              'qc_warnings': []},
                             {'position': 9864, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17,
                              'alt_count': 5, 'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
                              'qc_warnings': []},
                             {'position': 1895, 'ref_base': 'C', 'alt_base': 'A', 'ref_count': 17,
                              'alt_count': 5, 'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
                              'qc_warnings': []}]}}
    flanking = 15
    ref_file = os.path.join(os.getcwd(), "tests", "examples", "ref.fasta")
    for seq_record in SeqIO.parse(ref_file, "fasta"):
        ref_seq = "{}".format(seq_record.seq)
    case = __main__.add_tiles(scheme, ref_seq, flanking, conflict_pos)
    compare = {9: {}, 8: {1: [
        {'position': 359, 'ref_base': 'T', 'alt_base': 'C', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
         'qc_warnings': ['One or more degenerate bases'],
         'positive_tile': 'TGTTGAAGTTCGGCGCTACATCNGTGGCAAA',
         'negative_tile': 'TGTTGAAGTTCGGCGTTACATCNGTGGCAAA'},
        {'position': 515, 'ref_base': 'T', 'alt_base': 'G', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
         'qc_warnings': ['One or more degenerate bases'],
         'positive_tile': 'TCGGCGGNCAGGATGGTTTGCCGAATATCAG',
         'negative_tile': 'TCGGCGGNCAGGATGTTTTGCCGAATATCAG'},
        {'position': 9797, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': [],
         'positive_tile': 'CGGCAAAAATTTGCGTGATACCGCCGTAGAA',
         'negative_tile': 'CGGCAAAAATTTGCGCGATACCGCCGTAGAA'},
        {'position': 9864, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': [],
         'positive_tile': 'AAAACCGGCATTGTGTAGGTTAAGCAGAATG',
         'negative_tile': 'AAAACCGGCATTGTGCAGGTTAAGCAGAATG'},
        {'position': 1895, 'ref_base': 'C', 'alt_base': 'A', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8,
         'qc_warnings': ['One or more degenerate bases'],
         'positive_tile': 'GCCTGAATCTGGANAACTGGCAGGCGGAACT',
         'negative_tile': 'GCCTGAATCTGGANACCTGGCAGGCGGAACT'}]}}

    assert case == compare

