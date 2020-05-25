def tsv_to_membership(infile):
    groups = dict()
    with open(infile) as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        for row in reader:
            leaf = row[0]
            for i in range(1, len(row)):
                if row[i] == "":
                    if i== len(row)-1:
                        break
                    else:
                        print("There is an error in your tsv file")
                else:
                    if leaf not in groups.keys():
                        groups[leaf] = [row[i]]
                    else:
                        groups[leaf].append(row[i])

    return groups
