import csv
from argparse import ArgumentParser


def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser_inner = ArgumentParser(
        description='BIO_Hansel codes to group ')
    parser_inner.add_argument('--infile', type=str, required=True, help='file with sample names and biohansel codes')
    parser_inner.add_argument('--outfile', type=str, required=True, help='name for output file')
    return parser_inner.parse_args()


def main():
    args = parse_args()
    infile = args.infile
    outfile = args.outfile
    initial = []
    groups = dict()
    with open(infile) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        for row in reader:
            print(len(row))
            if len(row) ==0:
                continue
            if len(row) != 2:
                print("There is something wrong with the input file, "
                      "either a sample name is missing or a code is missing")
                raise SystemExit(0)
            leaf = row[0]
            bh_code = row[1]
            initial.append((leaf, bh_code))
    for i in range(0, len(initial)):
        to_split = initial[i][1]
        categories = to_split.split(".")
        last = ""
        for j in range(0,len(categories)):
            if j == 0:
                last = categories[0]
                add = last
                groups[initial[i][0]] = [add]
            else:
                add = last + "." + categories[j]
                last = add
                groups[initial[i][0]].append(add)
    out = open(f"{outfile}", "w+")
    for key in groups:
        line = f"{key}"
        for i in range(0,len(groups[key])):
            line = line + "\t" + groups[key][i]
        out.write(f"{line}\n")


# call main function
if __name__ == '__main__':
    main()