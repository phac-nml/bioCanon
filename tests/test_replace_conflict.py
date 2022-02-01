from bioCanon import __main__


def test_replace_conflict():
    conflict = 366
    start = 343
    scheme = {9: {}, 8: {1: [
        {'position': 359, 'ref_base': 'T', 'alt_base': 'C', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': [],
         'positive_tile': 'TGTTGAAGTTCGGCGCTACATCAGTGGCAAA',
         'negative_tile': 'TGTTGAAGTTCGGCGTTACATCAGTGGCAAA'},
        {'position': 515, 'ref_base': 'T', 'alt_base': 'G', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 9797, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 9864, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 1895, 'ref_base': 'C', 'alt_base': 'A', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []}]}}
    scheme_pieces = [8, 1, 0]
    __main__.replace_conflict(conflict, start, scheme, scheme_pieces)
    compare = {9: {}, 8: {1: [
        {'position': 359, 'ref_base': 'T', 'alt_base': 'C', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': [],
         'positive_tile': 'TGTTGAAGTTCGGCGCTACATCNGTGGCAAA',
         'negative_tile': 'TGTTGAAGTTCGGCGTTACATCNGTGGCAAA'},
        {'position': 515, 'ref_base': 'T', 'alt_base': 'G', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 9797, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 9864, 'ref_base': 'C', 'alt_base': 'T', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []},
        {'position': 1895, 'ref_base': 'C', 'alt_base': 'A', 'ref_count': 17, 'alt_count': 5,
         'positive_group': 'alt', 'g_id': 1, 'rank_id': 8, 'qc_warnings': []}]}}

    assert scheme == compare
