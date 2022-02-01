from bioCanon import __main__
import os
from ete3 import Tree


def test_parse_tree():
    tree_file = os.path.join(os.getcwd(),"tests", "examples", "testing.nwk")
    compare = Tree("((((((((((((Q,R),(O,P)),J),I),H),G),(S,F)),B),(C,(T,U))),"
                   "(((M,N),(K,L)),E)),(D,V)),A);")
    case = __main__.parse_tree(tree_file)
    assert print(case) == print(compare)

