from bioCanon import __main__
import os
import vcfpy


def test_alt_or_ref():
    vcf_file = os.path.join(os.getcwd(), "tests", "examples", "testing_single.vcf")
    reader = vcfpy.Reader.from_path(vcf_file)
    samples = reader.header.samples.names
    for record in reader:
        case = __main__.alt_or_ref(record, samples)
        compare = (['A', 'B', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q',
                    'R', 'S', 'V'], ['C', 'T', 'U'])
        assert case == compare
