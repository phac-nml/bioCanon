
"""
This tool is intended to wrap Biohansel and use it to provide summary Q&A for a
scheme developed in BioCannon.  It uses the same sorts of o
"""
import subprocess
import pandas as pd
from argparse import ArgumentParser
from collections import defaultdict
import os


def input_parser():
    """
    Parse the input arguments, use '-h' for help
        Returns
        -------
            arguments()
    """
    # Instantiate parser
    parser = ArgumentParser(description='BioQA - Scheme QA Tool')

    # Add arguments
    parser.add_argument('--in_scheme', type=str, required=True,
                        help='Expects a Biohansel scheme in .fasta format.')
    parser.add_argument('--in_seq_data', type=str, required=True,
                        help='Expects a directory containing .fasta assemblies or .fastq reads.')
    parser.add_argument('--biohansel_arguments', type=str, required=False, default='',
                        help='Pass as a quoted list any additional arguments to biohansel.')
    parser.add_argument('--filter_out_below', type=float, required=False, default=5.0,
                        help='Filter out kmers below the specified threshold.')
    parser.add_argument('--filter_out_above', type=float, required=False, default=1.0000001,
                        help='Filter out kmers above the specified threshold.')
    parser.add_argument('--output_kmer_name', type=str, required=False, default='bioqa_kmer.tsv',
                        help='Name of the file to store the QA kmer report.')
    parser.add_argument('--output_filtered_name', type=str, required=False, default='bioqa_filtered.tsv',
                        help='Name of the file to store a list of kmers which failed to meet the filtering criteria.')
    parser.add_argument('--outdir', type=str, required=False, default=os.getcwd(),
                        help='Output Directory to put results')
    # Return the parsed arguments
    return parser.parse_args()


if __name__ == '__main__':
    """ 
    Entry point for the script
    """
    # Set up a few filenames
    biohansel_kmer_fn = 'biohansel_output_kmer.txt'
    biohansel_summary_fn = 'biohansel_output_summary.txt'

    # Get the arguments
    inputs = input_parser()

    # Verify inputs files are present
    if not os.path.exists(inputs.in_scheme):
        print(f"{inputs.in_scheme} not found.  Input scheme required.")
        exit()

    if not os.path.exists(inputs.in_seq_data):
        print(f"{inputs.in_seq_data} not found.  Input sequence data required.")
        exit()

    # Verify that an acceptable version of Biohansel is installed
    cp = subprocess.run(['hansel', '-V'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    if cp.stderr:
        print(
            "There was an error attempting to run Biohansel (command 'hansel -V'). Please confirm that Biohanel is installed and you have permission to run it.")
        exit()

    # Run biohansel on the inputs
    cp = subprocess.run(
        ['hansel', '-s', inputs.in_scheme, '-D', inputs.in_seq_data, '--output-kmer-results', biohansel_kmer_fn,
         '--output-summary', biohansel_summary_fn, '--force', inputs.biohansel_arguments], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, shell=False)
    if cp.stderr:
        print("There was an error attempting to run Biohansel with the specified parameters.")
        exit()

    # Read in the input scheme
    # Note schemes are expected with labels in the format >[negative]numerical_position-Level_0.Level_i.Level_n
    pdict = defaultdict(dict)
    inf = open(inputs.in_scheme, 'r')
    line_cnt = 0
    for line in inf:
        if line[0] == '>':
            try:
                pos_line, hier_line = line[1:].split('-')
            except ValueError:
                print(f"Invalid label format (incorrect hyphenation) in line {line_cnt}.")
                continue

            # Parse the line
            if pos_line[:8] == 'negative':
                position = int(pos_line[8:])
                pdict[position]['negative'] = 0
                pdict[position]['nhierarchy'] = hier_line
            else:
                position = int(pos_line)
                pdict[position]['positive'] = 0
                pdict[position]['phierarchy'] = hier_line
        line_cnt += 1
    inf.close()

    if line_cnt % 2 != 0:
        print(f'Error! Uneven number of k-mers found.  Either a positive or negative tile is missing.')
    pos_cnt = line_cnt / 2
    print(f'{pos_cnt} k-mers described by the scheme.')

    # Check that each position has a positive and negative tile
    for pkey in pdict:
        if not ('positive' in pdict[pkey] and 'negative' in pdict[pkey]):
            print(f"Scheme Error! Position {pkey} does not have both a positive and negative tile.")
            continue
        if pdict[pkey]['phierarchy'] != pdict[pkey]['nhierarchy']:
            print(f"Scheme Error! Position {pkey} has inconsistent hierarchies for the positive and negative tiles.")

    # Check that the outputs exist
    if not os.path.exists(biohansel_kmer_fn):
        print(f"kmer summary '{biohansel_kmer_fn}' not generated and/or not accessible.")
        exit()

    if not os.path.exists(biohansel_summary_fn):
        print(f"summary '{biohansel_summary_fn}' not generated and/or not accessible.")
        exit()

    # Read in the kmer summary
    dfk = pd.read_csv(biohansel_kmer_fn, sep='\t')
    dfk_values = dfk['kmername'].value_counts()
    for el in dfk_values.index:
        pos_line = el.split('-')[0]
        if pos_line[:8] == 'negative':
            position = int(pos_line[8:])
            if position not in pdict:
                print(
                    f"Invalid tile in {biohansel_kmer_fn}: {position} not recognized as a position from {inputs.in_scheme}.")
                continue
            pdict[position]['negative'] += 1
        else:
            position = int(pos_line)
            if position not in pdict:
                print(
                    f"Invalid tile in {biohansel_kmer_fn}: {position} not recognized as a position from {inputs.in_scheme}.")
                continue
            pdict[position]['positive'] += 1

            # Read in the summary
    dfs = pd.read_csv(biohansel_summary_fn, sep='\t')
    # Get the number of samples
    num_samples = len(dfs)

    # Write the kmer report
    outf = open(os.path.join(os.path.join(os.getcwd(), inputs.outdir), inputs.output_kmer_name), "w")
    failed_pos = []
    outf.write(f"position\tpos tiles\tneg tiles\tmissing\t%present\n")
    for el in sorted(pdict.keys()):
        missing = num_samples - pdict[el]['positive'] - pdict[el]['negative']
        percent_present = float(missing) / float(num_samples)
        outf.write(f"{el}\t{pdict[el]['positive']}\t{pdict[el]['negative']}\t{missing}\t{percent_present}\n")
        if percent_present < inputs.filter_out_below or percent_present > inputs.filter_out_above:
            failed_pos.append(el)
    outf.close()

    # Write the set of kmers to exclude
    outf = open(os.path.join(os.path.join(os.getcwd(), inputs.outdir), inputs.output_filtered_name), "w")
    for el in failed_pos:
        outf.write(f"{el}\n")
    outf.close()

