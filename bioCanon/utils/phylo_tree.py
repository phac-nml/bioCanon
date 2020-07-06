from ete3 import Tree, NodeStyle, TreeStyle, TextFace,RectFace
import os



def parse_tree(tree_file,logging,ete_format=0,set_root=False,resolve_polytomy=True,ladderize=True,method='midpoint',outgroup=''):
    """
    Parses newick formatted tree into an ETE tree object
    Parameters
    ----------
    tree_file [str] : Path to newick formatted tree
    logging [logging obj] : logging object
    ete_format [int] : Refer to documentation http://etetoolkit.org/docs/latest/reference/reference_tree.html
    set_root [Bool] : Set root for the tree
    resolve_polytomy [Bool]: Force bifurcations of the tree on polytomies
    ladderize [Bool] : Sort the branches of a given tree (swapping children nodes) according to the size
    method [str]: Method to root tree, either midpoint or outgroup
    outgroup [str] : Name of taxon to root tree on

    Returns
    -------

    ETE tree obj

    """

    logging.info("Attempting to parse tree file {}".format(tree_file))

    if not os.path.isfile(tree_file):
        logging.error("Specified tree file {} is not found, please check that it exists".format(tree_file))
        return dict()
    if os.path.getsize(tree_file) == 0:
        logging.error("Specified tree file {} is found but is empty".format(tree_file))
        return dict()


    # Load a tree structure from a newick file.
    t = Tree(tree_file,format=ete_format)

    logging.info(
        "Read {} samples comprising {} nodes from tree {}".format((len(t.children[0].children[0].get_tree_root())),
                                                                  len(t.get_descendants()), tree_file))

    #Need this otherwise groups derived from the tree are inaccurate
    if resolve_polytomy:
        logging.info("Resolving polytomies through forced bifurcation {}".format(tree_file))
        t.resolve_polytomy()


    #Ladderize tree for aesthetics
    if ladderize:
        logging.info("Ladderizing tree {}".format(tree_file))
        t.ladderize()


    #Rooting tree based on user criteria
    if set_root:
        if method == 'midpoint':
            logging.info("Setting root based on midpoint rooting from ETE3")
            root = t.get_midpoint_outgroup()
            t.set_outgroup(root)
        elif method == 'outgroup':
            if outgroup == '':
                logging.error("User selected outgroup rooting but did not provide an outgroup, please refer to the documentation for setting an outgroup")
                return None
            else:
                #if outgroup label not in the tree, return an error to the user
                if t.search_nodes(name=outgroup):
                    t.set_outgroup(outgroup)
                else:
                    logging.error(
                        "Specified outgroup not found in tree, please check the name and try again")
                    return None

    sample_list = t.get_leaf_names()

    #label all internal nodes
    node_id = 0
    for node in t.traverse("preorder"):
        if node.name == '':
            while node_id in sample_list:
                node_id+=1
            node.name = str(node_id)
            node_id+=1


    num_samples = len(t.children[0].get_tree_root())
    num_nodes = len(t.get_descendants())
    is_rooted = len(t.get_tree_root().children) == 2
    root_id = t.get_tree_root().name
    max_node, max_dist = t.get_farthest_leaf()
    logging.info(
        "Read {} samples comprising {} nodes from tree {}:\n".format(num_samples,num_nodes, tree_file,sample_list))
    logging.info("Most distant sample id {}".format(max_node.name))
    logging.info("Max Distance {}".format(max_dist))
    if is_rooted:
        for node in t.traverse("levelorder"):
            if node.is_leaf():
                root_id = node.name
                break
        logging.info("Tree is rooted on sample {}".format(root_id))
    else:
        logging.error("Error the tree is unrooted, you must either specify the root using an outgroup or midpoint rooting")

    return t


def get_internal_clades(ete_tree_obj):
    """
    Identify internal nodes in the ETE3 tree and return a dict of the samples associated with that clade
    Parameters
    ----------
    ete_tree_obj [ETE tree obj] :

    Returns
    -------
    clade_dict [dict]: Dictionary of tree nodes which are not leaves

    """
    cache = ete_tree_obj.get_cached_content()
    clade_dict = {}
    for clade in cache:
        clade_id = clade.name
        children = clade.get_leaf_names()
        if len(children) == 1:
            continue
        clade_dict[clade_id] = children
    return clade_dict

def init_clade_info(ete_tree_obj):
    clade_info = {}
    memberships = get_internal_clades(ete_tree_obj)

    for clade_id in memberships:
        clade_info[clade_id] = {'num_members':len(memberships[clade_id]),'membership':set(memberships[clade_id]),'canSNPs':{},'canSNP_count':0}

    return clade_info

def visualize_tree(ete_tree_obj,genotypes={},metadata={},cansnp_supported_nodes={}):
    """
    Function for visualization of cannonical snp supported clades
    Parameters
    ----------
    ete_tree_obj
    cansnp_supported_nodes

    Returns
    -------

    """

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = True

    # Draws nodes as small red spheres of diameter equal to 10 pixels


    for n in ete_tree_obj.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "sphere"
        nstyle["fgcolor"] = "darkred"

        if n.name in cansnp_supported_nodes:
            print("Found")
            #nstyle["fgcolor"] = "blue"
            #nstyle["size"] = cansnp_supported_nodes[n.name]['size']
            #nstyle["bgcolor"] =cansnp_supported_nodes[n.name]['bgcolor']
            face = RectFace(width=100,height=100,fgcolor='red',bgcolor='green')


            n.add_face(face,column=1)
        n.set_style(nstyle)



    ete_tree_obj.show(tree_style=ts)