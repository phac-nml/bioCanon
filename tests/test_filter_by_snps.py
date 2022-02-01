from bioCanon import __main__


def test_filter_by_snps():
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
    case = __main__.filter_by_snps(scheme, 2)
    case.sort()
    compare = [359, 515, 1895, 9797, 9864]
    assert case == compare
