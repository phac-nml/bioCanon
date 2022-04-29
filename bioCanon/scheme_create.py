import os,shutil,sys
from argparse import ArgumentParser
import pandas as pd
sys.setrecursionlimit(1500)
from bioCanon.version import __version__
from bioCanon.utils import init_console_logger, read_fasta_dict, write_scheme
from bioCanon.utils.phylo_tree import parse_tree, get_tree_node_distances
from bioCanon.utils.vcf_helper import sanitize_vcf_file,pandas_vcf
from scipy.stats import entropy,fisher_exact
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

def parse_args():
    "Parse the input arguments, use '-h' for help"
    parser = ArgumentParser(
        description='bioCanon: Genotyping Scheme Creation Tool')

    parser.add_argument('--in_file', type=str, required=True, help='Either Variant Call SNP data (.vcf) or TSV SNP data (.txt)')
    parser.add_argument('--in_nwk', type=str, required=False, help='Newick Tree of strains')
    parser.add_argument('--in_groups', type=str, required=False, help='Tab delimited file of (sample_id|genotype)',
                        default=None)
    parser.add_argument('--in_meta', type=str, required=False, help='Tab delimited file of metadata', default=None)
    parser.add_argument('--reference', type=str, required=True, help='Reference fasta sequence from VCF')
    parser.add_argument('--outdir', type=str, required=True, help='Output Directory to put results')
    parser.add_argument('--prefix', type=str, required=False, help='Prefix for output files',default='bioCanon')
    parser.add_argument('--root_name', type=str, required=False, help='Name of sample to root tree',default='')
    parser.add_argument('--root_method', type=str, required=False, help='Method to root tree (midpoint,outgroup)', default=None)
    parser.add_argument('--bp_right', type=int, required=False, help='Num bases right of SNP to use for k-mer creation', default=8)
    parser.add_argument('--bp_left', type=int, required=False, help='Num bases downstream of SNP to use for k-mer creation',
                        default=8)
    parser.add_argument('--min_members', type=int, required=False, help='Minimum number of members for a clade to be valid',
                        default=1)
    parser.add_argument('--min_snps', type=int, required=False, help='Minimum number of unique snps for a clade to be valid',
                        default=1)
    parser.add_argument('--max_missing_position', type=float, required=False, help='Maximum percentage of samples which can be missing a position for it to be valid',
                        default=0.25)
    parser.add_argument('--max_missing_sample', type=float, required=False, help='Maximum percentage of positions which a sample can be missing for it to be included',
                        default=0.25)
    parser.add_argument('--min_polymorphism', type=float, required=False, help='Minimum number of samples with different states for a position to be valid',
                        default=0.0001)
    parser.add_argument('--fconst', type=str, required=False, help='Number of constant sites suitable for IQ-TREE -fconst',
                        default=None)
    parser.add_argument('--num_threads', type=int, required=False,
                        help='Number of threads to use',
                        default=1)
    parser.add_argument('--debug', required=False, help='Show debug information', action='store_true')
    parser.add_argument('-V', '--version', action='version', version="%(prog)s " + __version__)

    return parser.parse_args()

def calc_shanon_entropy(value_list):
    total = sum(value_list)
    if total == 0:
        return -1
    values = []
    for v in value_list:
        values.append(v / total)
    return entropy(values)

def calc_AMI(category_1, category_2):
    return adjusted_mutual_info_score(category_1, category_2, average_method='arithmetic')

def calc_ARI(category_1, category_2):
    return adjusted_rand_score(category_1, category_2)

def validate_args(arg_parse_obj,logging):
    """
    Provides validation of input parameters
    :param arg_parse_obj: [Argparse obj]
    :param logging: [logging ob]
    :return: [Bool] True if successful
    """
    success = True
    snp_data = arg_parse_obj.in_file

    if not os.path.isfile(snp_data):
        success = False
        logging.error("Error SNP data file {} not found".format(snp_data))
    elif os.path.getsize(snp_data) == 0:
        success = False
        logging.error("Error SNP data file {} is empty".format(snp_data))

    tree_file = arg_parse_obj.in_nwk

    if tree_file is not None:
        if not os.path.isfile(tree_file):
            success = False
            logging.error("Error tree file {} not found".format(tree_file))
        elif os.path.getsize(tree_file) == 0:
            success = False
            logging.error("Error tree file {} is empty".format(tree_file))

    group_file = arg_parse_obj.in_groups

    if group_file is not None:
        if not os.path.isfile(group_file):
            success = False
            logging.error("Error group file {} not found".format(tree_file))
        elif os.path.getsize(group_file) == 0:
            success = False
            logging.error("Error group file {} is empty".format(group_file))

    if (tree_file is None and group_file is None) or (tree_file is not None and group_file is not None):
        success = False
        logging.error("Error you need to specify either a tree file OR a group file ".format(tree_file))

    reference_fasta = arg_parse_obj.reference
    if not os.path.isfile(reference_fasta):
        success = False
        logging.error("Error fasta file {} not found".format(reference_fasta))
    elif os.path.getsize(reference_fasta) == 0:
        success = False
        logging.error("Error fasta file {} is empty".format(reference_fasta))

    out_dir = arg_parse_obj.outdir


    prefix = arg_parse_obj.prefix
    invalid_chars = [' ','/',"\\"]
    for char in prefix:
        if str(char) in invalid_chars:
            success = False
            logging.error("Error prefix contains invalid characters [{}]".format(prefix,",".join(invalid_chars)))

    method = arg_parse_obj.root_method
    outgroup = arg_parse_obj.root_name
    valid_methods = ['midpoint','outgroup']
    if method not in valid_methods:
        success = False
        logging.error("Error bioCanon only supports midpoint or outgroup for tree rooting, you specified: {}".format(method))
    elif method == 'outgroup' and outgroup == '':
        success = False
        logging.error("Error you selected outgroup method for rooting but left --root_name blank, please specify the taxon name to root the tree")

    return success

def identify_canonical_snps(snp_data,ete_tree_obj):
    """
    Examines each SNP and determines which clade in the tree exclusively contains all of the members positive for the
    SNP (if possible)
    :param snp_data: [dict]
    :param ete_tree_obj: [ETE obj]
    :return: snp_data [dict] with clades that are canonical
    """
    snp_clade_info = {

    }
    tree_samples = ete_tree_obj.get_leaf_names()

    for chrom in snp_data:
        if not chrom in snp_clade_info:
            snp_clade_info[chrom] = {}
        for position in snp_data[chrom]:
            snp_clade_info[chrom][position] = {}
            ambig_samples = {}
            if 'N' in snp_data[chrom][position]:
                ambig_samples = ambig_samples.update(snp_data[chrom][position]['N'])
            if '-'in snp_data[chrom][position]:
                ambig_samples = ambig_samples.update(snp_data[chrom][position]['-'])

            for base in snp_data[chrom][position]:
                if base == 'N' or base == '-':
                    continue
                snp_clade_info[chrom][position][base] = {
                    'clade_id': '',
                    'clade_total_members': 0,
                    'clade_ambig_members': 0,
                    'clade_postitive_members':0,
                    'clade_negative_members':0,
                    'is_canonical': False,
                    'out_bases': []
                }

                snp_samples = snp_data[chrom][position][base]['samples']


                if len(snp_samples) == 1:
                    ancestor_node = ete_tree_obj.get_leaves_by_name(snp_samples[0])[0].get_ancestors()[0]
                    node_leaves = snp_samples
                else:
                    ancestor_node = ete_tree_obj.get_common_ancestor(snp_samples)
                    node_leaves = ancestor_node.get_leaf_names()

                node_leaves_noambig = set(node_leaves) - set(ambig_samples)

                snp_clade_info[chrom][position][base]['clade_id'] = ancestor_node.name
                snp_clade_info[chrom][position][base]['clade_total_members'] = len(node_leaves)
                snp_clade_info[chrom][position][base]['clade_ambig_members'] = len(node_leaves) - len(node_leaves_noambig)
                snp_clade_info[chrom][position][base]['clade_postitive_members'] = len(snp_samples)
                snp_clade_info[chrom][position][base]['clade_negative_members'] = len(list(set(node_leaves_noambig) - set(snp_samples)))
                snp_clade_info[chrom][position][base]['out_bases'] = list( set(snp_data[chrom][position][base]['out_bases']) - set('N'))
                if set(node_leaves) == set(snp_samples) or set(snp_samples) == set(node_leaves_noambig):
                    snp_clade_info[chrom][position][base]['is_canonical'] = True
                else:
                    if len(snp_data[chrom][position][base]['out_bases']) == 1:
                        snp_samples = list(set(tree_samples) - set(snp_samples))
                        ancestor_node = ete_tree_obj.get_common_ancestor(snp_samples)
                        node_leaves = ancestor_node.get_leaf_names()
                        node_leaves_noambig = set(node_leaves) - set(ambig_samples)
                        if set(node_leaves) == set(snp_samples) or set(snp_samples) == set(node_leaves_noambig):
                            del(snp_clade_info[chrom][position][base])
                            out_bases = [base]
                            base = snp_data[chrom][position][base]['out_bases'][0]
                            snp_clade_info[chrom][position][base] = {
                                'clade_id': ancestor_node.name,
                                'clade_total_members': len(node_leaves),
                                'clade_ambig_members': len(node_leaves) - len(node_leaves_noambig),
                                'clade_postitive_members': len(snp_samples),
                                'clade_negative_members': len(list(set(node_leaves_noambig) - set(snp_samples))),
                                'is_canonical': True,
                                'out_bases': out_bases
                            }

    return snp_clade_info

def summarize_node_support(snp_clade_info):
    """
    Conversts SNP data into a clade indexed datastructure
    :param snp_clade_info: [dict] Dictionary of clades and SNPs supporting them
    :return: [dict] Summary of clade SNP support
    """
    clades = {}
    for chrom in snp_clade_info:
        for position in snp_clade_info[chrom]:
            for base in snp_clade_info[chrom][position]:
                data = snp_clade_info[chrom][position][base]
                clade_id = data['clade_id']
                if clade_id not in clades:
                    clades[clade_id] = {'count_canonical':0,'count_homoplasy':0,'clade_total_members':data['clade_total_members'],'canonical_snps':{}}
                if data['is_canonical']:
                    clades[clade_id]['count_canonical']+= 1
                    if not chrom in clades[clade_id]['canonical_snps']:
                        clades[clade_id]['canonical_snps'][chrom] = {}
                    clades[clade_id]['canonical_snps'][chrom][position] = {'in_base':base,'out_bases':data['out_bases']}
                else:
                    clades[clade_id]['count_homoplasy'] += 1
    return clades

def remove_unsupported_clades(ete_tree_obj, valid_clades):
    """
    Remove non-leaf nodes from the tree while maintaining children
    :param ete_tree_obj: [ETE obj] Phylogenetic tree
    :param valid_clades: [list] Valid node names to maintain
    :return: ete_tree_obj: [ETE obj] Phylogenetic tree with only nodes in the valid list
    """
    tree = ete_tree_obj.copy()

    for node in tree.traverse("preorder"):
        if node.is_root() or node.is_leaf():
            continue

        node_id = node.name
        remove_node = True
        if node_id in valid_clades:
            remove_node = False
        if remove_node:
            node.delete()

    return tree

def get_count_valid_internal_nodes(node_children,clade_snp_support,min_members=1,min_snps=1):
    """
    Checks node children to see if they are internal tree nodes supported by canonical SNPs
    :param node_children: [ETE obj] Of nodes
    :param clade_snp_support: [dict] Dictionary indexed by nodes in the tree with SNP information
    :param min_members: [int] Minimum number of clade members for clade to be considered valid
    :param min_snps: [int] Minimum number of canonical SNPs supporting clade
    :return:
    """
    count_internal_nodes = 0
    for n in node_children:

        if n.is_leaf() or n.name not in clade_snp_support:
            continue

        count_cansnp = clade_snp_support[n.name]['count_canonical']
        count_members = clade_snp_support[n.name]['clade_total_members']

        if count_cansnp < min_snps or count_members < min_members:
            continue

        count_internal_nodes += 1
    return count_internal_nodes

def get_valid_nodes(ete_tree_obj,clade_snp_support,min_members=1,min_snps=1):
    """

    :param ete_tree_obj: [ETE obj] Tree object
    :param clade_snp_support: [dict] Dictionary indexed by nodes in the tree with SNP information
    :param min_members: [int] Minimum number of clade members for clade to be considered valid
    :param min_snps: [int] Minimum number of canonical SNPs supporting clade
    :return: [List] of nodes supported by cannonical SNPs and meeting group size
    """
    tree = ete_tree_obj
    valid_nodes = []

    for node in tree.traverse("levelorder"):
        node_id = node.name

        if node.is_leaf() :
            continue
        node_children = node.get_children()
        count_internal_nodes = get_count_valid_internal_nodes(node_children,clade_snp_support,min_members,min_snps)

        if node_id in clade_snp_support:
            clade_snp_support[node_id]['count_valid_child_nodes'] = count_internal_nodes
            count_cansnp = clade_snp_support[node_id]['count_canonical']
            count_members = clade_snp_support[node_id]['clade_total_members']

            if count_members < min_members or count_cansnp < min_snps:
                print("{} {} {} skip".format(node_id,count_members, count_cansnp ))
                continue
        else:
            continue
        valid_nodes.append(node_id)
    print(sorted(valid_nodes))
    return valid_nodes

def prune_tree(ete_tree_obj,valid_nodes):
    """
    :param ete_tree_obj: [ETE obj] Tree object
    :param valid_nodes: [List] Node names which are valide
    :return: [ETE Tree ] With only valid nodes
    """
    invalid_nodes = []
    for node in ete_tree_obj.traverse("postorder"):
        node_id = node.name
        if node.is_leaf() or node_id not in valid_nodes:
            continue
        node_ancestors = node.get_ancestors()

        for n in node_ancestors:
            ancestor_id = n.name
            if n.is_leaf() or \
                    ancestor_id not in valid_nodes or \
                    ancestor_id in invalid_nodes:
                continue
            children = n.get_children()
            count_internal_nodes = 0
            for c in children:
                if c.is_leaf() or c.is_root() or c.name in invalid_nodes:
                    continue
                count_internal_nodes += 1
            if count_internal_nodes < 1:
                invalid_nodes.append(c.name)
    valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
    pruned_tree = remove_unsupported_clades(ete_tree_obj, valid_nodes)
    return [pruned_tree,valid_nodes]

def get_internal_nodes(genotypes):
    terminal_nodes = count_terminal_nodes(genotypes)
    internal_nodes = []
    for sample_id in genotypes:
        geno = genotypes[sample_id]
        length = len(geno)
        if length == 1:
            continue
        for i in range(1,length):
            node = geno[i]
            if not node in terminal_nodes and not node in internal_nodes:
                internal_nodes.append(node)
    return internal_nodes

def get_invalid_terminal_nodes(terminal_nodes,threshold=1):
    invalid_nodes = []
    for t_node in terminal_nodes:
        if terminal_nodes[t_node] < threshold:
            invalid_nodes.append(t_node)
    return invalid_nodes

def count_terminal_nodes(genotypes):
    terminal_nodes = {}
    for sample_id in genotypes:
        geno = genotypes[sample_id]
        #Never remove the root
        if len(geno) <= 1:
            continue
        t_node = geno[len(geno)-1]
        if not t_node in terminal_nodes:
            terminal_nodes[t_node]=0
        terminal_nodes[t_node] += 1
    return terminal_nodes

def generate_genotypes(tree):
    samples = tree.get_leaf_names()
    geno = {}
    for sample_id in samples:
        geno[sample_id] = []
    for node in tree.traverse("preorder"):
        node_id = node.name
        if node.is_leaf() or node.is_root():
            continue
        leaves = node.get_leaf_names()
        for sample_id in leaves:
            geno[sample_id].append(node_id)
    return geno

def compress_genotypes(ete_tree_obj,clade_snp_support,min_members=1,min_snps=1,min_perc=0.0,compress_geno=False):
    valid_nodes = get_valid_nodes(ete_tree_obj, clade_snp_support, min_members, min_snps)
    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, valid_nodes)
    genotypes = generate_genotypes(pruned_tree)
    if not compress_geno:
        return genotypes

    invalid_nodes = get_invalid_terminal_nodes(count_terminal_nodes(genotypes), min_members)
    valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
    pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)

    # get rid of terminal nodes which don't meet minimum group size
    while (len(invalid_nodes) > 0):
        genotypes = generate_genotypes(pruned_tree)
        invalid_nodes = get_invalid_terminal_nodes(count_terminal_nodes(genotypes), min_members)
        valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
        pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)

    genotypes = generate_genotypes(pruned_tree)
    print(genotypes)
    division_info = identify_division_entropies(genotypes)
    entropies = division_info['entropies']
    num_samples = len(genotypes)
    sample_counts = division_info['sample_counts']
    # filter 0 entropy divisions, do not remove the "root" node
    # filter ranks where the perc of samples being subdivided is less than threshold
    invalid_ranks = []
    percs = []
    for i in range(1, len(entropies)):
        e = entropies[i]
        if e == 0:
            invalid_ranks.append(i)
        p = sample_counts[i] / num_samples
        if p < min_perc:
            invalid_ranks.append(i)
        percs.append(p)

    invalid_nodes = set()
    for sample_id in genotypes:
        for r, node in enumerate(genotypes[sample_id]):
            if r in invalid_ranks:
                invalid_nodes.add(genotypes[sample_id][r])

    valid_nodes = list(set(valid_nodes) - invalid_nodes)
    pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)
    genotypes = generate_genotypes(pruned_tree)
    print(genotypes)
    #filter out individual nodes which do not meet the min_perc
    node_counts = {}
    for sample_id in genotypes:
        for node in genotypes[sample_id]:
            if not node in node_counts:
                node_counts[node] =0
            node_counts[node] +=1

    invalid_nodes = set()
    for node in node_counts:
        p = node_counts[node] / num_samples
        if p < min_perc:
            invalid_nodes.add(node)

    valid_nodes = list(set(valid_nodes) - invalid_nodes)
    pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)

    genotypes = generate_genotypes(pruned_tree)
    invalid_nodes = get_internal_nodes(genotypes)

    # get rid of internal only nodes
    while (len(invalid_nodes) > 0):
        genotypes = generate_genotypes(pruned_tree)
        invalid_nodes = get_internal_nodes(genotypes)
        valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
        pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)

    return generate_genotypes(pruned_tree)

def identify_division_entropies(genotype_assignments):
    max_level = 0

    #get max depth of genotypes
    for sample_id in genotype_assignments:
        num_levels = len(genotype_assignments[sample_id])
        if num_levels > max_level:
            max_level = num_levels

    entropies = []
    sample_counts = []
    for i in range(0, max_level):
        sample_count = 0
        genotypes = {}
        rank = i + 1
        for sample_id in genotype_assignments:
            #skip genotypes which are not resolved at this level
            if len(genotype_assignments[sample_id]) < rank:
                continue
            genotype = '.'.join(genotype_assignments[sample_id][0:rank])
            if not genotype in genotypes:
                genotypes[genotype] = 0
            genotypes[genotype]+=1
            sample_count += 1
        sample_counts.append(sample_count)
        entropies.append(calc_shanon_entropy(list(genotypes.values())))

    return {'entropies':entropies,'sample_counts':sample_counts}

def associate_genotype_snps(sample_genotypes,clade_snp_support,snp_clade_info):
    genotypes = []
    genotype_snps = {}

    for sample_id in sample_genotypes:

        geno = sample_genotypes[sample_id]
        #print("{}\t{}".format(sample_id, ".".join(geno)))

        if geno not in genotypes:
            genotypes.append(geno)
        else:
            continue

        for clade_id in sample_genotypes[sample_id]:
            if clade_id not in clade_snp_support:
                continue

            cansnps = clade_snp_support[clade_id]['canonical_snps']

            for chrom in cansnps:

                for pos in cansnps[chrom]:
                    if not chrom in genotype_snps:
                        genotype_snps[chrom] = {}

                    in_base = cansnps[chrom][pos]['in_base']
                    out_bases = cansnps[chrom][pos]['out_bases']

                    for i in range(0,len(sample_genotypes[sample_id])):

                        if sample_genotypes[sample_id][i] == clade_id:
                            geno = ".".join(sample_genotypes[sample_id][:i+1])
                            genotype_snps[chrom][pos] = {'in_base':in_base,'out_bases':out_bases,'genotype':geno,'clade_id':clade_id}

    return genotype_snps

def tree_based_scheme_generation(vcf_file,tree_file,metadata, prefix,outdir,logging,root_method,outgroup='',min_members=1, min_snps=1,max_alt_states=4, disruptive_threshold=1,max_missing=0.25,compress_geno=False):
    """
    Generate a bioHansel compatible SNP typing scheme based on a phylogenetic tree
    :param vcf_file:
    :param tree_file:
    :param prefix:
    :param outdir:
    :param logging:
    :param root_method:
    :param outgroup:
    :param min_members:
    :param min_snps:
    :param max_alt_states:
    :param disruptive_threshold:
    :return:
    """

    #Report files
    cleaned_vcf_file = os.path.join(outdir,"{}.cleaned.vcf".format(prefix))
    snp_report_file = os.path.join(outdir, "{}.original.snp.report.txt".format(prefix))
    sample_report_file = os.path.join(outdir, "{}.original.sample.report.txt".format(prefix))


    ete_tree_obj =  parse_tree(tree_file,logging,ete_format=1,set_root=True,resolve_polytomy=True,ladderize=True,method=root_method,outgroup=outgroup)
    tree_samples = ete_tree_obj.get_leaf_names()
    (snp_samples,snp_data) = sanitize_vcf_file(vcf_file, cleaned_vcf_file, snp_report_file, sample_report_file, logging, min_members, max_missing,
                      max_alt_states=max_alt_states, disruptive_threshold=disruptive_threshold)

    #Check that the samples in both files are the same and if not then report the error
    intersection = list(set(tree_samples) & set(snp_samples))

    if len(intersection) != len(tree_samples) or len(intersection) != len(snp_samples):
        missing = list(set(snp_samples) - set(tree_samples))
        if len(missing) > 0:
            logging.warning("Your SNP data file {} contains the following samples which are not in your tree file {}:\n {}".format(vcf_file,tree_file,
                list(set(snp_samples) - set(tree_samples))
            ))

        missing = list(set(tree_samples) - set(snp_samples))
        if len(missing) > 0:
            logging.error("Your tree data file {} contains the following samples which are not in your snp file {}:\n {}".format(tree_file,vcf_file,
                list(set(tree_samples) - set(snp_samples))
            ))

    snp_clade_info = identify_canonical_snps(snp_data, ete_tree_obj)
    canonical_snp_positions = []


    #get canonical positions
    for chrom in snp_clade_info:
        for pos in snp_clade_info[chrom]:
            for base in snp_clade_info[chrom][pos]:
                if snp_clade_info[chrom][pos][base]['is_canonical']:
                    canonical_snp_positions.append(pos)
                    break

    cansnp_info = {}
    canonical_node_ids = []
    for chrom in snp_clade_info:
        for pos in snp_clade_info[chrom]:
            cansnp_info[pos] = {
                'chrom':chrom,
                'pos':pos,
                'A_count':0,
                'T_count': 0,
                'C_count': 0,
                'G_count': 0,
                'N_count': 0,
                'N_count': 0,
                'A_node_id': '',
                'T_node_id': '',
                'C_node_id': '',
                'G_node_id': '',
                'N_node_id': '',
                'N_node_id': '',
            }

            for base in snp_clade_info[chrom][pos]:
                count_key = "{}_count".format(base)
                node_key = "{}_node_id".format(base)
                clade_id = snp_clade_info[chrom][pos][base]['clade_id']
                canonical_node_ids.append(clade_id)
                clade_total_members = snp_clade_info[chrom][pos][base]['clade_total_members']
                cansnp_info[pos][count_key] = clade_total_members
                cansnp_info[pos][node_key] = clade_id
    cansnp_outfile = os.path.join(outdir,"{}-cansnp.info.txt".format(prefix))
    pd.DataFrame.from_dict(cansnp_info,orient='index').to_csv(cansnp_outfile,sep="\t",header=True,index=False)

    #Write output file showing the original tree representation of the samples
    genotypes = generate_genotypes(ete_tree_obj)
    original_tree_genotype_outfile = os.path.join(outdir, "{}-original-tree-genotypes.txt".format(prefix))
    fh = open(original_tree_genotype_outfile,'w')
    fh.write("sample_id\tgenotype\n")
    for sample_id in genotypes:
        genotype = ".".join(genotypes[sample_id])
        fh.write("{}\t{}\n".format(sample_id,genotype))
    fh.close()

    clade_snp_support = summarize_node_support(snp_clade_info)
    sample_genotypes = compress_genotypes(ete_tree_obj, clade_snp_support, min_members, min_snps,compress_geno)


        #Write output file showing the modified tree representation of the samples
    mod_tree_genotype_outfile = os.path.join(outdir, "{}-mod-tree-genotypes.txt".format(prefix))
    fh = open(mod_tree_genotype_outfile,'w')
    fh.write("sample_id\tgenotype\n")
    for sample_id in sample_genotypes:
        genotype = ".".join(sample_genotypes[sample_id])
        fh.write("{}\t{}\n".format(sample_id,genotype))
    fh.close()

    associations = metadata_association(metadata, sample_genotypes)

    return associate_genotype_snps(sample_genotypes,clade_snp_support,snp_clade_info)

def sliding_window_search(reference_sequence,pos,kmer_right_bp=8,kmer_left_bp=8,min_length=17):
    kmer_right_bp = 8
    kmer_left_bp = 8
    min_length = 17
    ref_len = len(reference_sequence)
    target_length = kmer_left_bp + kmer_right_bp + 1
    nt_count = 0
    ambig_count = 0
    end = pos


    start = pos - 1

    #find initial start
    for i in range(0,ref_len):
        base = reference_sequence[start].upper()
        if base in ['A', 'T', 'C', 'G']:
            nt_count += 1
        if base in ['N','Y','W','R','V','B']:
            ambig_count += 1

        if start < 0:
            start = 0
            break
        if nt_count == min_length:
            if base not in ['A', 'T', 'C', 'G']:
                start -=1
            break
        start -= 1
        if nt_count >= target_length:
            break

    prev_start = -1
    prev_end = pos
    prev_window_ambig_count = 99999999

    for i in range(0,min_length):
        nt_count = 0
        ambig_count = 0

        for k in range(start,end):
            base = reference_sequence[k].upper()
            if base in ['A', 'T', 'C', 'G']:
                nt_count += 1
            if base in ['N', 'Y', 'W', 'R', 'V', 'B']:
                ambig_count += 1

        if ambig_count < prev_window_ambig_count and nt_count >= min_length:
            if reference_sequence[end] == 'N' and end != pos:
                end -=1
            if reference_sequence[start] == 'N':
                start -=1
            prev_start = start
            prev_window_nt_count = nt_count
            prev_window_ambig_count = ambig_count
        start -= 1
        if (nt_count + ambig_count) >= target_length:
            break

    start = prev_start

    for i in range(0,(pos-start)):
        nt_count = 0
        ambig_count = 0

        for k in range(start,end):
            base = reference_sequence[k].upper()
            if base in ['A', 'T', 'C', 'G']:
                nt_count += 1
            if base in ['N', 'Y', 'W', 'R', 'V', 'B']:
                ambig_count += 1

        if ambig_count < prev_window_ambig_count and nt_count >= min_length:
            if reference_sequence[end] == 'N' and end != pos:
                end -=1
            if reference_sequence[start] == 'N':
                start -=1
            prev_start = start
            prev_end = end
            prev_window_nt_count = nt_count
            prev_window_ambig_count = ambig_count
        start += 1
        end += 1
        if (nt_count + ambig_count) >= target_length:
            break

    if start < 0:
        prev_start = 0
    if end < pos :
        prev_end = pos

    return (prev_start,prev_end)

def count_kmers(seq, K=2):
    mers = {}
    """Count kmers in sequence"""
    for i in range(0,len(seq)):
        mer = list(seq[i:i+K])
        mer.sort()
        mer = ''.join(mer)
        if len(mer) != K:
            continue
        if not mer in mers:
            mers[mer] = 0
        mers[mer]+=1
    return mers

def find_initial_start(pos,reference_sequence,min_length):
    """
    Using an initial position, it finds the initial starting position which satisfies the minimum length
    :param pos: int
    :param reference_sequence: str
    :param min_length: int
    :return: int
    """
    ref_len=len(reference_sequence)
    nt_count = 0
    start = pos - 1
    for i in range(0,ref_len):
        base = reference_sequence[start].upper()
        if base in ['A', 'T', 'C', 'G']:
            nt_count += 1
        if start < 0:
            start = 0
            break
        if nt_count >= min_length:
            break
        start -= 1
    return start

def optimize_kmer(pos,reference_sequence,min_length,max_length,max_ambig=5,min_complexity=0.2):
    """
    Accepts a position and a sequence and determines the kmer stretch which maximizes length, complexity and minimizes
    ambiguous characters
    :param pos: int position
    :param reference_sequence: str reference sequence
    :param min_length: int minimum length of kmer
    :param max_length: int maximum length of kmer
    :param max_ambig:  int maximum number of iupac characters
    :param min_complexity: float maximum percentage composition of one 2-mer
    :return:
    """

    prev_score = 0
    opt_kmer = [-1,-1]
    rlen = len(reference_sequence)
    start = find_initial_start(pos, reference_sequence, max_length) +1

    for length_target in range(min_length,max_length):

        for k in range(start ,pos):

            s = pos - k
            if s > length_target :
                continue
            valid = True
            rel_start = k
            nt_count = 0
            rel_end = k
            base_count = 0
            while nt_count < length_target:
                if base_count >= length_target or rel_end >= rlen -1:
                    break
                base = reference_sequence[rel_end]
                if base in ['A', 'T', 'C', 'G']:
                    nt_count += 1
                if base != '-':
                    base_count+=1
                rel_end += 1

            if start <= 0 or start >= rlen or rel_end >= rlen or rel_end < pos:
                continue
            kmer = reference_sequence[rel_start:rel_end].replace('-','')
            klen = len(kmer)


            if klen > max_length :
                continue

            #count ambiguous characters
            bases = ['A','T','C','G']
            nt_count = 0
            mers = count_kmers(kmer, K=1)
            for b in bases:
                if b in mers:
                    nt_count+=mers[b]
            count_ambig = klen - nt_count
            if count_ambig > max_ambig:
                valid = False

            #determine the complexity of the sequence and remove kmers composed heavily of the same 2-mer
            mers = count_kmers(kmer, K=2)
            num_mers = sum(mers.values())
            mer_perc = []

            for m in mers:
                if mers[m]/num_mers > min_complexity:
                    valid = False
                    break
                mer_perc.append(mers[m]/num_mers )

            if not valid:
                continue

            minimum = min(mer_perc)
            if max_ambig > 0:
                score = 5*(1 - ((nt_count+count_ambig)/max_length)) + (1 - minimum) + (1 - count_ambig/max_ambig)
            else:
                score = 5*(1 - ((nt_count+count_ambig)/max_length)) + (1 - minimum) + 1

            if prev_score < score:
                prev_score = score
                opt_kmer = [rel_start,rel_end]

    return opt_kmer

def generate_kmers(reference_sequence,genotype_snps,kmer_right_bp=12,kmer_left_bp=12,min_len=15,max_len=21,max_ambig=4):
    """
    Accepts a fasta reference sequence and then extracts out the sequence window of interest for the genotyping kmer
    Then format the kmer into in or out k-mers in a format that can be printed later
    :param reference_sequence: [dict] Fasta dict of sequences indexed by id
    :param genotype_snps: [dict] Fasta dict of genotyping SNPs
    :param kmer_right_bp: [int] Select this many bp to the right of the SNP position
    :param kmer_left_bp: [int] Select this many bp to the left of the SNP position
    :return: [dict] of Sequence kmers
    """
    kmers = {}
    iupac_lookup ={
        'A':'A',
        'T':'T',
        'C':'C',
        'G':'G',
        'AT':'W',
        'AC':'M',
        'AG':'R',
        'CT':'Y',
        'GT':'K',
        'CG':'S',
        'CGT':'B',
        'AGT':'D',
        'ACT':'H',
        'ACG':'V'
    }
    negated_bases = {
        'A':'B',
        'T':'V',
        'C':'D',
        'G':'H'
    }

    start_positions = []
    for chrom in genotype_snps:
        if chrom not in reference_sequence:
            if len(reference_sequence) == 1:
                key = list(reference_sequence.keys())[0]
                ref_len = len(reference_sequence[key])
                ref_seq = generate_consensus(reference_sequence[key], genotype_snps[chrom])
            else:
                break
        else:
            ref_len = len(reference_sequence[chrom])
            ref_seq = generate_consensus(reference_sequence[chrom], genotype_snps[chrom])
        if not chrom in kmers:
            kmers[chrom] = {}
        sorted_positions = sorted(list(genotype_snps[chrom].keys()))
        for pos in sorted_positions:
            if not pos in kmers[chrom]:
                kmers[chrom][pos] = {'genotype':'','in_kmer':'','out_kmers':[],'kmer_start':0,'kmer_end':0}

            genotype = genotype_snps[chrom][pos]['genotype']
            in_base = genotype_snps[chrom][pos]['in_base']
            out_bases =  genotype_snps[chrom][pos]['out_bases']

            kmers[chrom][pos]['genotype'] = genotype

            (start, end) = optimize_kmer(pos, ref_seq, min_len, max_len, max_ambig, min_complexity=0.6)
            #fall back to other approach if no suitable ones were found
            if start == -1 or end == -1:
                (start,end) = sliding_window_search(ref_seq, pos)

            if start == -1 or end == -1:
                continue

            while start in start_positions and start < pos:
                start+= 1
                end+=1
                if end >= ref_len:
                    end = ref_len-1
            start_positions.append(start)



            kmers[chrom][pos]['kmer_start'] = start
            kmers[chrom][pos]['kmer_end'] = end

            seq = list(ref_seq[start:end])

            #Modify the SNP position base

            snp_pos =  pos - start -1
            if len(seq) < snp_pos:
                continue
            seq[snp_pos] = in_base
            in_kmer = ''.join(seq)
            kmers[chrom][pos]['in_kmer'] = "".join(seq).replace('-','')
            out_bases.sort()
            iupac = ''.join(out_bases)
            if iupac in iupac_lookup:
                seq[snp_pos] = iupac_lookup[iupac]
            else:
                seq[snp_pos] = negated_bases[in_base]
            out_kmer = ''.join(seq)
            seq = ref_seq[start:end]
            kmers[chrom][pos]['out_kmers'].append("".join(seq).replace('-',''))
    return kmers

def metadata_association(metadata,genotype_assignments):
    if len(metadata) == 0:
        return {}
    fields = list(metadata[list(metadata.keys())[0]].keys())
    genotypes = sorted(list(set(genotype_assignments.values())))
    max_depth = 0
    for genotype in genotypes:
        value = len(genotype.split('.'))
        if value > max_depth:
            max_depth = value
    for sample_id in genotype_assignments:
        genotype_assignments[sample_id] = genotype_assignments[sample_id] .split('.')
    field_values = {}
    field_value_counts = {}
    genotype_counts = {}
    total_samples = len(genotype_assignments)
    associations = {}
    for sample_id in genotype_assignments:
        for i in range(0, max_depth):
            if not i in associations:
                associations[i] = {}
            genotype = '.'.join(genotype_assignments[sample_id][0:i + 1])
            if not i in genotype_counts:
                genotype_counts[i] = {}
                field_value_counts[i] = {}
            if not genotype in genotype_counts[i]:
                genotype_counts[i][genotype] = {'total': 0, 'fields': {}}
            for field in fields:
                value = metadata[sample_id][field]
                if not field in associations[i]:
                    associations[i][field] = {'entropy': {}, 'ari': {}, 'ami': {}, 'fischer': {}}
                if not field in genotype_counts[i][genotype]['fields']:
                    genotype_counts[i][genotype]['fields'][field] = {}
                if not value in genotype_counts[i][genotype]['fields'][field]:
                    genotype_counts[i][genotype]['fields'][field][value] = 0
                genotype_counts[i][genotype]['fields'][field][value]+=1

                if field not in field_value_counts[i]:
                    field_value_counts[i][field] = {}
                if not value in field_value_counts[i][field]:
                    field_value_counts[i][field][value] = 0
                field_value_counts[i][field][value]+=1
            genotype_counts[i][genotype]['total'] += 1

    for i in range(0,max_depth):
        field_values[i] = {}
        for field in fields:
            field_value_counts[field] = {}
            field_values[i][field] = {}
            category_1 = []
            category_2 = []
            for sample_id in genotype_assignments:
                value = metadata[sample_id][field]
                if not value in field_value_counts[field]:
                    field_value_counts[field][value] = 0
                field_value_counts[field][value]+=1
                genotype = '.'.join(genotype_assignments[sample_id][0:i + 1])
                category_1.append(value)
                category_2.append(genotype)
                if not genotype in field_values[i][field]:
                    field_values[i][field][genotype] = {}
                if not value in field_values[i][field][genotype]:
                    field_values[i][field][genotype][value] = 0
                field_values[i][field][genotype][value] += 1
            associations[i][field]['entropy'] = calc_shanon_entropy(list(field_values[i][field][genotype].values()))
            associations[i][field]['ari'] = calc_ARI(category_1, category_2)
            associations[i][field]['ami'] = calc_AMI(category_1, category_2)

    for i in range(0, max_depth):
        for genotype in genotype_counts[i]:
            genotype_total = genotype_counts[i][genotype]['total']
            for field in genotype_counts[i][genotype]['fields']:
                if not genotype in associations[i][field]['fischer']:
                    associations[i][field]['fischer'][genotype] = {}

                for value in genotype_counts[i][genotype]['fields'][field]:

                    pos_geno_count = genotype_counts[i][genotype]['fields'][field][value]
                    neg_geno_count = genotype_total - genotype_counts[i][genotype]['fields'][field][value]
                    pos_field_count = field_value_counts[i][field][value] - pos_geno_count
                    neg_field_count = total_samples - (pos_geno_count + neg_geno_count + pos_field_count)

                    table = [
                        [pos_geno_count, neg_geno_count],
                        [pos_field_count,neg_field_count]
                    ]
                    oddsr, p = fisher_exact(table, alternative='greater')
                    associations[i][field]['fischer'][genotype][value] = {'odds_ratio':oddsr,'p':p}

    return associations

def generate_consensus(seq,genotype_snps):
    iupac_lookup ={
        'A':'A',
        'T':'T',
        'C':'C',
        'G':'G',
        'AT':'W',
        'AC':'M',
        'AG':'R',
        'CT':'Y',
        'GT':'K',
        'CG':'S',
        'CGT':'B',
        'AGT':'D',
        'ACT':'H',
        'ACG':'V'
    }
    negated_bases = {
        'A':'B',
        'T':'V',
        'C':'D',
        'G':'H'
    }
    seq = list(seq)
    for pos in genotype_snps:
        in_base = genotype_snps[pos]['in_base']
        out_bases = genotype_snps[pos]['out_bases']
        bases = list(set(list(in_base) + out_bases + list(seq[pos])))
        bases = list(set(bases) - set(['N']))

        bases.sort()
        key = ''.join(bases)
        if key in iupac_lookup:
            char = iupac_lookup[key]
        else:
            char = 'N'
        seq[pos] = char

    return ''.join(seq)

def parse_metadata(file):
    df = pd.read_csv(file,sep="\t",header=0)
    columns = set(df.columns.tolist()) - set('sample_id')
    metadata = {}
    for index, row in df.iterrows():
        metadata[row['sample_id']] = {}
        for field in columns:
            metadata[row['sample_id']][field] = row[field]
    return  metadata


def run():
    cmd_args = parse_args()
    logging = init_console_logger(3)
    status = validate_args(cmd_args,logging)
    if not status:
        logging.error("Cannot continue, please read error messages and corret them before continuing")
        sys.exit()

    #Inputs
    snp_data_file = cmd_args.in_file
    in_meta = cmd_args.in_meta
    outdir = cmd_args.outdir
    tree_file = cmd_args.in_nwk
    root_method = cmd_args.root_method
    outgroup = cmd_args.root_name
    prefix = cmd_args.prefix
    reference = cmd_args.reference
    bp_right = cmd_args.bp_right
    bp_left = cmd_args.bp_left
    min_members = cmd_args.min_members
    min_snps = cmd_args.min_snps

    if in_meta is not None:
        metadata = parse_metadata(in_meta)
    else:
        metadata = {}

    if not os.path.isdir(outdir):
        logging.info("Creating analysis directory {}".format(outdir))
        os.mkdir(outdir, 0o755)

    #Derrived filenames
    genotyping_snps = tree_based_scheme_generation(snp_data_file, tree_file, metadata, prefix, outdir, logging, root_method, outgroup,min_members,min_snps)
    reference_seq_dict = read_fasta_dict(reference)

    biohansel_scheme_outfile = os.path.join(outdir, "{}-biohansel.scheme.fasta".format(prefix))
    logging.info("Writting bioHansel compatible scheme to {}".format(biohansel_scheme_outfile))
    write_scheme( generate_kmers(reference_seq_dict, genotyping_snps, kmer_right_bp=bp_right, kmer_left_bp=bp_left),biohansel_scheme_outfile)
