from bioCanon import __main__
import os


def test_tree_path():
    group_info = "none"
    tree_file = os.path.join(os.getcwd(), "tests", "examples", "testing.nwk")
    case = __main__.path_check(group_info, tree_file)
    compare = ({'Q': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                'R': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2],
                'O': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3],
                'P': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4],
                'J': [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0, 0],
                'I': [1, 1, 1, 1, 1, 1, 1, 1, 2, 0, 0, 0],
                'H': [1, 1, 1, 1, 1, 1, 1, 2, 0, 0, 0, 0],
                'G': [1, 1, 1, 1, 1, 1, 2, 0, 0, 0, 0, 0],
                'S': [1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0],
                'F': [1, 1, 1, 1, 1, 2, 4, 0, 0, 0, 0, 0],
                'B': [1, 1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0],
                'C': [1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0],
                'T': [1, 1, 1, 2, 4, 3, 0, 0, 0, 0, 0, 0],
                'U': [1, 1, 1, 2, 4, 4, 0, 0, 0, 0, 0, 0],
                'M': [1, 1, 2, 3, 5, 5, 0, 0, 0, 0, 0, 0],
                'N': [1, 1, 2, 3, 5, 6, 0, 0, 0, 0, 0, 0],
                'K': [1, 1, 2, 3, 6, 7, 0, 0, 0, 0, 0, 0],
                'L': [1, 1, 2, 3, 6, 8, 0, 0, 0, 0, 0, 0],
                'E': [1, 1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0],
                'D': [1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'V': [1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'A': [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]},
               ['A', 'D', 'V', 'E', 'B', 'C', 'T', 'U', 'M', 'N', 'K', 'L', 'G', 'S', 'F', 'H',
                'I', 'J', 'Q', 'R', 'O', 'P'], 'tree')
    assert case == compare


def test_groups_path():
    tree_file = "none"
    group_info = os.path.join(os.getcwd(), "tests", "examples", "testing.tsv")
    case = __main__.path_check(group_info, tree_file)
    compare = ({'A': ['1', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'B': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.1', 0, 0, 0, 0, 0, 0, 0],
                'C': ['2', '2.1', '2.1.1', '2.1.1.2', '2.1.1.2.1', 0, 0, 0, 0, 0, 0, 0],
                'D': ['2', '2.2', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                'E': ['2', '2.1', '2.1.2', '2.1.2.1', '2.1.1.1.2', 0, 0, 0, 0, 0, 0, 0],
                'F': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.2', 0, 0, 0, 0, 0, 0],
                'G': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.2',
                      0, 0, 0, 0, 0],
                'H': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.2', 0, 0, 0, 0],
                'I': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.2', 0, 0, 0],
                'J': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.1', '2.1.1.1.2.1.1.1.1.2', 0, 0],
                'K': ['2', '2.1', '2.1.2', '2.1.2.2', '2.1.2.2.2', 0, 0, 0, 0, 0, 0, 0],
                'L': ['2', '2.1', '2.1.2', '2.1.2.2', '2.1.2.2.2', 0, 0, 0, 0, 0, 0, 0],
                'M': ['2', '2.1', '2.1.2', '2.1.2.2', '2.1.2.2.1', 0, 0, 0, 0, 0, 0, 0],
                'N': ['2', '2.1', '2.1.2', '2.1.2.2', '2.1.2.2.1', 0, 0, 0, 0, 0, 0, 0],
                'O': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.1', '2.1.1.1.2.1.1.1.1.1',
                      '2.1.1.1.2.1.1.1.1.1.2', 0],
                'P': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.1', '2.1.1.1.2.1.1.1.1.1',
                      '2.1.1.1.2.1.1.1.1.1.2', 0],
                'Q': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1', '2.1.1.1.2.1.1',
                      '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.1', '2.1.1.1.2.1.1.1.1.1',
                      '2.1.1.1.2.1.1.1.1.1.1', 0],
                'R': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.1',
                      '2.1.1.1.2.1.1', '2.1.1.1.2.1.1.1', '2.1.1.1.2.1.1.1.1',
                      '2.1.1.1.2.1.1.1.1.1', '2.1.1.1.2.1.1.1.1.1.1', 0],
                'S': ['2', '2.1', '2.1.1', '2.1.1.1', '2.1.1.1.2', '2.1.1.1.2.2', 0, 0, 0, 0, 0, 0],
                'T': ['2', '2.1', '2.1.1', '2.1.1.2', '2.1.1.2.2', 0, 0, 0, 0, 0, 0, 0],
                'U': ['2', '2.1', '2.1.1', '2.1.1.2', '2.1.1.2.2', 0, 0, 0, 0, 0, 0, 0],
                'V': ['2', '2.2', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]},
               ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
                'Q', 'R', 'S', 'T', 'U', 'V'], 'groups')
    assert case == compare
