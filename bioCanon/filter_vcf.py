"""
During the generation of Biohansel schemes filtering of the vcf files may be necessary to remove
snp locations that are not helpful for a variety of reasons.  This script filters a vcf file to keep
 or remove desired positions and also outputs a fasta file with all the snp position values
 concatenated for the samples so that a new, refined phylogeny can be constructed
"""
import os
import csv
from argparse import ArgumentParser
import vcfpy


def parse_args():
    """
    Parse the input arguments, use '-h' for help

        Returns
        -------
            parser_inner.parse_args()
    """
    parser_inner = ArgumentParser(
        description='BioHansel Scheme development accessory script')
    parser_inner.add_argument('--in_vcf', type=str, required=True, help='vcf File of SNPs')
    parser_inner.add_argument('--outdir', type=str, required=False, default=os.getcwd(),
                              help='Output Directory to put results')
    parser_inner.add_argument('--keep_csv', type=str, required=False,
                              help='csv file with comma deliminated positions to be kept for the'
                                   ' reduced fasta file')
    parser_inner.add_argument('--remove_csv', type=str, required=False,
                              help='csv file with comma deliminated positions to be removed for the'
                                   ' reduced fasta file')
    parser_inner.add_argument('--o', type=str, required=False, default="Reduced",
                              help='name for output files')
    return parser_inner.parse_args()


def main():
    """
    The main body of the code that calls all the functions required to filter the vcf files and
    generate fasta files

    """
    args = parse_args()
    vcf_file = args.in_vcf
    out = args.o
    out_path = args.outdir
    try:
        reader = vcfpy.Reader.from_path(vcf_file)
    except FileNotFoundError:
        print(f"There was an error in reading {vcf_file}.  Verify that the file exists")
        raise SystemExit(0)
    except vcfpy.exceptions.IncorrectVCFFormat:
        print(f"There was an error in reading {vcf_file}.  Verify that the file is in the correct "
              f"format")
        raise SystemExit(0)
    if not args.keep_csv:
        infile = args.remove_csv
        opp = 0
    elif not args.remove_csv:
        infile = args.keep_csv
        opp = 1
    else:
        print("Select one reduction method: either list positions to keep or list positions to "
              "remove")
        raise SystemExit(0)
    try:
        pos_list = []
        with open(infile) as tsvfile:
            keep = csv.reader(tsvfile, delimiter=",")
            for record in keep:
                pos_list.append(record)
    except FileNotFoundError:
        print("There was an error reading the csv file.  Check the file name and try again")
        raise SystemExit(0)

    fasta_dict = dict()
    samples = reader.header.samples.names
    writer = vcfpy.Writer.from_path(f'{out}.vcf', reader.header)
    in_vcf = []
    for record in reader:
        # skip over the metadata lines at the top of the file
        if not record.is_snv():
            continue
        position = f"{record.POS}"
        in_vcf.append(position)
        if opp == 1:
            if position in pos_list:
                writer.write_record(record)
                # pull the positional information and the reference value and all the alt values
                # from the vcf file
                line = [record.CHROM, record.POS, record.REF]
                line += [alt.value for alt in record.ALT]

                # initialize variables and pull the ref and alt SNPs
                alt_base = record.ALT[0].value
                ref_base = record.REF
                tracker = 0

                # go through the samples and divide them by alt vs ref bases
                for call in record.calls:
                    state = int(call.data.get('GT'))
                    sample_name = samples[tracker]
                    if sample_name not in fasta_dict.keys():
                        fasta_dict[sample_name] = ""
                    if state == 0:
                        fasta_dict[sample_name] += ref_base
                    else:
                        fasta_dict[sample_name] += alt_base
                    tracker += 1
        if opp == 0:
            if position not in pos_list:
                writer.write_record(record)
                # pull the positional information and the reference value and all the alt values
                # from the vcf file
                line = [record.CHROM, record.POS, record.REF]
                line += [alt.value for alt in record.ALT]

                # initialize variables and pull the ref and alt SNPs
                alt_base = record.ALT[0].value
                ref_base = record.REF
                tracker = 0

                # go through the samples and divide them by alt vs ref bases
                for call in record.calls:
                    state = int(call.data.get('GT'))
                    sample_name = samples[tracker]
                    if sample_name not in fasta_dict.keys():
                        fasta_dict[sample_name] = ""
                    if state == 0:
                        fasta_dict[sample_name] += ref_base
                    else:
                        fasta_dict[sample_name] += alt_base
                    tracker += 1
    # display a warning if the supplied csv list is not a subset of all the positions in the file
    if not all(x in in_vcf for x in pos_list):
        print("WARNING:\nOne or more positions in the supplied csv file are not found in the "
              "original vcf file")
    flag = 0
    in_vcf.sort()
    for i in range(0, (len(in_vcf)-1)):
        if in_vcf[i] == in_vcf[i+1]:
            flag = 1

    if flag:
        print("WARNING:\nOne or more entries in the supplied vcf file contained the same numerical"
              " position value.  Both entries were processed in the desired manner i.e. if the "
              "position value was to be kept both were kept, if they were to be removed they were "
              "both removed")

    fasta_file = open(os.path.join(out_path, f"{out}.fasta"), "w+")
    for key in fasta_dict.keys():
        fasta_file.write(f">{key}\n")
        fasta_file.write(f"{fasta_dict[key]}\n")


if __name__ == '__main__':
    main()
