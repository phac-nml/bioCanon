import os,shutil,sys
from argparse import ArgumentParser
from bioCanon.version import __version__
from bioCanon.utils import init_console_logger, read_fasta_dict
from bioCanon.utils.phylo_tree import parse_tree
from bioCanon.utils.vcf_helper import sanitize_vcf_file

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
    #if os.path.isfile(out_dir) or os.path.isdir(out_dir):
    #    success = False
    #    logging.error("Error analysis directory {} exists, please specify a new location".format(out_dir))
    #else:
     #   # initialize analysis directory
     #   if not os.path.isdir(out_dir):
      #      os.mkdir(out_dir, 0o755)

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

    for chrom in snp_data:
        if not chrom in snp_clade_info:
            snp_clade_info[chrom] = {}
        for position in snp_data[chrom]:
            snp_clade_info[chrom][position] = {}
            if 'N' in snp_data[chrom][position]:
                ambig_samples = snp_data[chrom][position]['N']
            else:
                ambig_samples = []

            for base in snp_data[chrom][position]:
                if base == 'N':
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
                ancestor_node = ete_tree_obj.get_common_ancestor(snp_samples)

                node_leaves = ancestor_node.get_leaf_names()
                node_leaves_noambig = set(node_leaves) - set(ambig_samples)

                snp_clade_info[chrom][position][base]['clade_id'] = ancestor_node.name
                snp_clade_info[chrom][position][base]['clade_total_members'] = len(node_leaves)
                snp_clade_info[chrom][position][base]['clade_ambig_members'] = len(node_leaves) - len(node_leaves_noambig)
                snp_clade_info[chrom][position][base]['clade_postitive_members'] = len(snp_samples)
                snp_clade_info[chrom][position][base]['clade_negative_members'] = len(list(set(node_leaves_noambig) - set(snp_samples)))
                snp_clade_info[chrom][position][base]['out_bases'] = snp_data[chrom][position][base]['out_bases']
                if set(node_leaves) == set(snp_samples) or set(node_leaves) == set(node_leaves_noambig):
                    snp_clade_info[chrom][position][base]['is_canonical'] = True



    return snp_clade_info

def summarize_node_support(snp_clade_info):
    print(snp_clade_info)
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
                continue
        else:
            continue

        # Do not assign a code to an internal node which does not split samples
        if count_internal_nodes == 0:
            continue
        valid_nodes.append(node_id)

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
        if len(geno) == 1:
            continue
        t_node = geno[len(geno)-1]
        if not t_node in terminal_nodes:
            terminal_nodes[t_node]=0
        terminal_nodes[t_node] += 1
    return terminal_nodes

def generate_genotypes(tree):
    geno = {}
    for node in tree.traverse("preorder"):
        node_id = node.name
        if node.is_leaf() or node.is_root():
            continue
        leaves = node.get_leaf_names()
        for sample_id in leaves:
            if sample_id not in geno:
                geno[sample_id] = []
            geno[sample_id].append(node_id)
    return geno

def compress_genotypes(ete_tree_obj,clade_snp_support,min_members=1,min_snps=1):
    valid_nodes = get_valid_nodes(ete_tree_obj, clade_snp_support, min_members, min_snps)
    pruned_tree, valid_nodes = prune_tree(ete_tree_obj, valid_nodes)

    genotypes = generate_genotypes(pruned_tree)
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
    invalid_nodes = get_internal_nodes(genotypes)

    # get rid of internal only nodes
    while (len(invalid_nodes) > 0):
        genotypes = generate_genotypes(pruned_tree)
        invalid_nodes = get_internal_nodes(genotypes)
        valid_nodes = list(set(valid_nodes) - set(invalid_nodes))
        pruned_tree, valid_nodes = prune_tree(pruned_tree, valid_nodes)

    return generate_genotypes(pruned_tree)

def associate_genotype_snps(sample_genotypes,clade_snp_support,snp_clade_info):
    genotypes = []
    genotype_snps = {}
    for sample_id in sample_genotypes:
        geno = ".".join(sample_genotypes[sample_id])
        if geno not in genotypes:
            genotypes.append(geno)
        else:
            continue

        for clade_id in sample_genotypes[sample_id]:
            if clade_id not in clade_snp_support:
                continue
            cansnps = clade_snp_support[clade_id]['canonical_snps']

            print(cansnps)
            for chrom in cansnps:
                print(cansnps[chrom])
                for pos in cansnps[chrom]:
                    if not chrom in genotype_snps:
                        genotype_snps[chrom] = {}

                    in_base = cansnps[chrom][pos]['in_base']
                    out_bases = cansnps[chrom][pos]['out_bases']

                    #Remove Ambiguous and in_base

                    genotype_snps[chrom][pos] = {'in_base':in_base,'out_bases':out_bases,'genotype':geno,'clade_id':clade_id}

    return genotype_snps










def tree_based_scheme_generation(vcf_file,tree_file,prefix,outdir,logging,root_method,outgroup='',min_members=1, min_snps=1,max_alt_states=1, disruptive_threshold=1,max_missing=0.25):
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
    node_report_file = os.path.join(outdir, "{}.original.nodes.report.txt".format(prefix))

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

    #print(snp_clade_info )

    clade_snp_support = summarize_node_support(snp_clade_info)
    #print(clade_snp_support)

    sample_genotypes = compress_genotypes(ete_tree_obj, clade_snp_support, min_members, min_snps)
    print(sample_genotypes)
    return associate_genotype_snps(sample_genotypes,clade_snp_support,snp_clade_info)






def generate_kmers(reference_sequence,genotype_snps,kmer_right_bp=15,kmer_left_bp=15):
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
    for chrom in genotype_snps:
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

            start = pos - kmer_right_bp - 1
            end = pos + kmer_left_bp

            kmers[chrom][pos]['kmer_start'] = start
            kmers[chrom][pos]['kmer_end'] = end

            kmer_bases = range(start,end)
            #print(reference_sequence.keys())
            seq = list(reference_sequence[chrom][start:end])

            # determine if there are other bases which need to be masked
            for i in kmer_bases:
                if i in sorted_positions:
                    index = i - start -1
                    seq[index] = 'N'

            #Modify the SNP position base
            snp_pos =  pos - start -1
            seq[snp_pos] = in_base
            kmers[chrom][pos]['in_kmer'] = "".join(seq)
            for b in out_bases:
                if b == 'N':
                    continue
                seq[snp_pos] = b
                kmers[chrom][pos]['out_kmers'].append("".join(seq))

    return kmers




def run():
    cmd_args = parse_args()
    logging = init_console_logger(3)
    status = validate_args(cmd_args,logging)
    if not status:
        logging.error("Cannot continue, please read error messages and corret them before continuing")
        sys.exit()

    #Inputs
    snp_data = cmd_args.in_file
    outdir = cmd_args.outdir
    tree_file = cmd_args.in_nwk
    root_method = cmd_args.root_method
    outgroup = cmd_args.root_name
    prefix = cmd_args.prefix
    reference = cmd_args.reference
    reference_seq_dict = read_fasta_dict(reference)

    genotyping_snps = tree_based_scheme_generation(snp_data, tree_file, prefix, outdir, logging, root_method, outgroup)



    print(generate_kmers(reference_seq_dict, genotyping_snps, kmer_right_bp=15, kmer_left_bp=15))