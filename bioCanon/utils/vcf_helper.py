import vcfpy, sys, os
import numpy as np
import pandas as pd

def get_vcf_row_calls(record, format_id, samples, logging):
    """

    Parameters
    ----------
    record [vcfpy reader record obj] : must be provided by an iterating vcpy reader obj
    format_id [str] : vcf format id for genotype
    samples [list] : list of samples in vcf file

    Returns
    -------

    data [dict] : indexed on position and containing the chromosome and sequences associated with each state in vcf

    """
    chromosome = record.CHROM
    position = record.POS
    ref_base = record.REF

    data = {position: {'bases': {}, 'chrom': chromosome, 'ref_base': ref_base}}

    # More extensive handelling required for multi-state alt-bases
    alt = record.ALT
    alt_bases = []
    count_alt_bases = 0
    if len(alt) > 1:
        for i in range(0, len(alt)):
            if alt[i].value in ('-', '*', 'X', 'x'):
                alt_bases.append('N')
            else:
                count_alt_bases += 1
                alt_bases.append(alt[i].value)
    else:
        alt_bases.append(alt[0].value)

    length = len(record.calls)

    # Process each sample and add it to the base it is associated with
    for i in range(0, length):
        call = record.calls[i]
        sample_id = samples[i]

        if not format_id in call.data:
            logging.error(
                "Format ID : {} not found in vcf file base calls {}, check file is valid vcf ".format(format_id,
                                                                                                      call.data))
            return []

        base = call.data[format_id]
        if isinstance(base, int) or base.isnumeric():
            base = int(base)
            if base == 0:
                base = ref_base
            else:
                base = alt_bases[base - 1]

        if not base in data[position]['bases']:
            data[position]['bases'][base] = []

        data[position]['bases'][base].append(sample_id)

    return data


def vcf_to_fasta(vcf_file, out_fasta, logging, position_filter={}, sample_filter=[]):
    """

    Parameters
    ----------
    vcf_file [str] : Path to vcf file
    out_fasta [str] : Path to write results to
    logging [logging obj]: Logging object
    position_filter [dict] : dictionary indexed by {chrom : [positions]} to exclude from fasta
    sample_filter [list] : list of sample_ids to exclude from fasta

    Returns
    -------
    Success [bool] : True if success
    """
    logging.info("Reading vcf file {}".format(vcf_file))
    reader = vcfpy.Reader.from_path(vcf_file)
    samples = reader.header.samples.names
    format_id = reader.header.format_ids()[0]
    sequences = {}

    # initialize set of samples
    for sample_id in samples:
        if sample_id in sample_filter:
            continue
        sequences[sample_id] = {}

    logging.info("Processing each vcf call".format(vcf_file))
    for record in reader:
        calls = get_vcf_row_calls(record, format_id, samples)
        position = calls.keys()[0]
        chrom = calls['chrom']
        bases = calls['bases']

        if chrom in position_filter and position in position_filter[chrom]:
            continue

        for base in bases:
            for sample_id in bases[base]:
                if not sample_id in sequences:
                    continue
                if not chrom in sequences[sample_id]:
                    sequences[sample_id][chrom] = []
                sequences[sample_id][chrom].append(base)

    if write_fasta_dict(sequences, out_fasta):
        logging.info("Writting sequences to {}...success".format(out_fasta))
        return True
    else:
        logging.error("Writting sequences to {}...failed".format(out_fasta))
        return False


def sanitize_vcf_file(vcf_file, out_file, snp_log_file, sample_log_file, logging, min_count=1, max_missing=0.25,
                      max_alt_states=1, disruptive_threshold=1,window_size=30,max_snps=2):
    """
    Filter a user provided vcf and write a filtered vcf file
    Parameters
    ----------
    vcf_file [str] : Path to vcf file
    out_file [str] : Path to write santized vcf file
    snp_log_file [str] : Path to write SNP log report
    sample_log_file [str]: Path to write sample report
    logging [logging obj] : logging object
    min_count [int] : Minimum number of variable samples ie. if set to 2, at least two samples must differ from the rest
    max_missing [float] : Percent of samples which can be missing data for a position for that position to be valid
    max_alt_states [int] : Maxmimum number of bases which can be in the vcf for that position to be valid
    disruptive_threshold [int]: Identify sequences which are blunting resolution by looking at which ones are missing
                                where other sequences have it. By default, it looks for sequences which are uniquely
                                missing a position conserved in all of the other samples

    Returns
    -------

    List of sample_id's in vcf

    """
    logging.info("Reading vcf file {}".format(vcf_file))
    reader = vcfpy.Reader.from_path(vcf_file)
    header = reader.header
    samples = header.samples.names
    num_samples = len(samples)
    snp_report = open(snp_log_file, 'w')
    sample_report = open(sample_log_file, 'w')

    sample_ambig_counts = {}
    for sample in samples:
        sample_ambig_counts[sample] = {'total_ambig': 0, 'uniq_ambig': 0, 'low_freq_variants': 0}
    format_id = header.format_ids()[0]

    logging.info("VCF file contains #{} samples".format(len(samples)))
    logging.info("VCF file samples: {}".format(samples))

    logging.info("Writing filtered vcf file to {}".format(out_file))
    writer = vcfpy.Writer.from_path(out_file, reader.header)

    count_postitions = 0
    valid_position_count = 0
    filtered_positions_count = 0
    filtered_positions_ambig = 0
    filtered_positions_alt = 0
    ambiguous_positions = {}
    snp_data = {}
    positions = []
    for record in reader:
        calls = get_vcf_row_calls(record, format_id, samples, logging)
        count_postitions += 1
        position = list(calls.keys())[0]
        positions.append(position)
        chromosome = calls[position]['chrom']
        bases = calls[position]['bases']
        ref_base = calls[position]['ref_base']

        if not chromosome in snp_data:
            snp_data[chromosome] = {}
        if not position in snp_data[chromosome]:
            snp_data[chromosome][position] = {}

        count_ref = 0
        count_alt = 0
        count_missing = 0
        count_alt_bases = 0
        alt_samples = []

        for base in bases:

            if base == ref_base:
                count_ref = len(bases[base])
            elif base == 'N':
                count_missing = len(bases[base])
                if not position in ambiguous_positions:
                    ambiguous_positions[position] = []
                ambiguous_positions[position] = bases[base]
                for sample_id in ambiguous_positions[position]:
                    sample_ambig_counts[sample_id]['total_ambig'] += 1
            else:
                count_alt_bases += 1
                count_alt += len(bases[base])
                alt_samples += bases[base]


        if count_ref >= min_count and count_alt >= min_count and count_missing / num_samples <= max_missing and count_alt_bases <= max_alt_states:
            snp_data[chromosome][position][base] = {'samples':bases[base],'out_bases':list(set(list(bases.keys())) - set([base]))}
            valid_position_count += 1
            writer.write_record(record)
        else:
            if count_ref < min_count or count_alt < min_count:
                filtered_positions_count += 1
                snp_report.write(
                    "Filtering due to minimum polymorphic sample count\tchrom: {}\tposition: {}\tcount_ref: {}\tcount_alt: {}\tcount_missing: {}\n".format(
                        chromosome, position, count_ref, count_alt, count_missing))
                for sample_id in alt_samples:
                    sample_ambig_counts[sample_id]['low_freq_variants'] += 1
            if count_missing / num_samples >= max_missing:
                filtered_positions_ambig += 1
                snp_report.write(
                    "Filtering due to maximum missing sample count\tchrom: {}\tposition: {}\tcount_ref: {}\tcount_alt: {}\tcount_missing: {}\n".format(
                        chromosome, position, count_ref, count_alt, count_missing))
            if count_alt_bases > max_alt_states:
                filtered_positions_alt += 1
                snp_report.write(
                    "Filtering due to more than allowed alt base states\tchrom: {}\tposition: {}\tcount_ref: {}\tcount_alt: {}\tcount_missing: {}\n".format(
                        chromosome, position, count_ref, count_alt, count_missing))

    #for position in positions:

    logging.info("Read {} positions in file {}".format(count_postitions, vcf_file))
    logging.info("Filtered {} positions due to minimum polymorphic sample requirement".format(filtered_positions_count))
    logging.info("Filtered {} positions due to missing in more than {}% of samples".format(filtered_positions_ambig,
                                                                                           max_missing))
    logging.info("Filtered {} positions due to more than allowed alternative base states".format(filtered_positions_alt,
                                                                                                 max_alt_states))
    logging.info("{} positions are valid after filtering".format(valid_position_count))

    # Identify for the user unusual sequences that they might want to remove from their analyses
    disruptive_sequence_check = identify_disruptive_sequences(ambiguous_positions, disruptive_threshold)

    for sample_id in disruptive_sequence_check:
        sample_ambig_counts[sample_id]['unique_ambig'] = disruptive_sequence_check[sample_id]

    # Get 95 percentile for each attribute
    percentile_total_ambig = get_percentile(sample_ambig_counts, 'total_ambig', percentile=95)
    percentile_unique_ambig = get_percentile(sample_ambig_counts, 'uniq_ambig', percentile=95)
    percentile_low_freq_var = get_percentile(sample_ambig_counts, 'low_freq_variants', percentile=95)

    logging.info(
        "95% Percentile: Total Ambig={}\t95% Percentile: Unique Ambig={}\t95% Percentile: low frequency variants={}".format(
            percentile_total_ambig, percentile_unique_ambig, percentile_low_freq_var))

    for sample_id in disruptive_sequence_check:
        status = 'PASS'
        if sample_ambig_counts[sample_id]['total_ambig'] / count_postitions > max_missing:
            status = 'FAIL'
        elif sample_ambig_counts[sample_id]['unique_ambig'] > percentile_unique_ambig:
            status = 'WARNING: Sample has unusually high number of unique missing core positions'
        elif sample_ambig_counts[sample_id]['low_freq_variants'] > percentile_low_freq_var:
            status = 'WARNING: Sample has unusually high number of low frequency variants'
        elif sample_ambig_counts[sample_id]['total_ambig'] > percentile_total_ambig:
            status = 'WARNING: Sample has unusually high number of missing core positions'

        sample_report.write("{}\tTOTAL_AMBIG={}\tUNIQUE_AMBIG={}\tLOW_FREQ_VARIANTS={}\tSTATUS={}\n".format(sample_id,
                                                                                                            sample_ambig_counts[
                                                                                                                sample_id][
                                                                                                                'total_ambig'],
                                                                                                            sample_ambig_counts[
                                                                                                                sample_id][
                                                                                                                'unique_ambig'],
                                                                                                            sample_ambig_counts[
                                                                                                                sample_id][
                                                                                                                'low_freq_variants'],
                                                                                                            status))

    sample_report.close()
    snp_report.close()
    return [samples, snp_data]

def identify_snp_positions_to_mask(positions,window_size=30,max_snps=2):
    """
    Uses a window approach to identify high density areas of SNPs which will generate highly degenrate k-mers
    :param positions: [List] Integer positions to check if there are too many
    :param window_size: [int] Window size to consider too many SNPs
    :param max_snps: [int] Mazimum allowable SNPs in a windo
    :return: [Lis] of positions to mask
    """
    mask = []
    for i in len(positions):
        pos = positions[i]
        count = 0
        window = [pos]
        for k in range(pos+1,pos+window_size):
            if k in positions:
                count+=1
                window.append(k)
        if count > max_snps:
            mask += window
    return mask





    return


def get_percentile(sample_dict, key, percentile=10):
    """
    Get the perctile for a numeric list of samples
    Parameters
    ----------
    sample_dict [dict] : {sample_id : {'key': [numeric values] }
    key [str] : Key of
    percentile [int] : 1 - 99 Percentile desired

    Returns
    -------
    percentile value

    """
    values = []
    for id in sample_dict:
        if key in sample_dict[id]:
            values.append(sample_dict[id][key])
    if len(values) == 0:
        return None
    return np.percentile(np.array(values), percentile)


def identify_disruptive_sequences(ambiguous_positions, threshold):
    """
    Count the number of times a given sample introduces a missing character into a position which would be otherwise core
    Parameters
    ----------
    ambiguous_positions [dict] : {position: [sample list]}
    threshold [int] : Max number of sequences required to be missing sequence for it to be considered disruptive

    Returns
    -------
    disruptive_samples[dict] : sample_id and number of times a sequence was disruptive
    """
    disruptive_samples = {}
    for pos in ambiguous_positions:
        total = len(ambiguous_positions[pos])
        if total < threshold:
            for sample_id in ambiguous_positions[pos]:
                if sample_id not in disruptive_samples:
                    disruptive_samples[sample_id] = 0
                disruptive_samples[sample_id] += 1
    return disruptive_samples


def write_fasta_dict(fasta_dict, out_file):
    """
    Write sequence dictionary in fasta format
    Parameters
    ----------
    fasta_dict [dict]: sample_id indexed dictionary of sequences
    out_file [str]: Path of file to write sequences to

    Returns
    -------

    Success [bool] : True if write was a success

    """

    try:
        fasta = open(out_file, 'w')
    except:
        return False
    for sample_id in fasta_dict:
        seq = []
        for seq_id in fasta_dict[sample_id]:
            seq += fasta_dict[sample_id][seq_id]
        fasta.write(">{}\n{}\n".format(sample_id, "".join(seq)))
    fasta.close()
    return True


class pandas_vcf:
    required_fields = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    def find_header_line(self, file):
        line_num = 0
        with open(file) as search:
            for line in search:
                if line[0:5] == '#CHROM':
                    return [line_num,line]
                line_num += 1
        return [-1,'']

    def find_missing_fields(self):
        header = self.get_header()
        missing_fields = []
        for field in self.required_fields:
            if field not in header:
                missing_fields.append(field)
        return missing_fields

    def set_header(self,header):
        self.header = header

    def get_header(self):
        return self.header

    def set_sample_list(self,samples):
        self.sample_list = samples

    def get_sample_list(self):
        return self.sample_list

    def check_file(self, file):
        if not os.path.isfile(file):
            return False
        if os.path.getsize(file) == 0:
            return False
        return True


    def __init__(self,file):
        if not self.check_file(file):
            return False
        (line_num,header_line) = self.find_header_line()
        self.set_header(header_line.split("\t"))
        missing_fields = self.find_missing_fields()
        if len(missing_fields) > 0:
            return False
        self.set_sample_list(list(set(self.get_header()) - set(self.required_fields)))
        if line_num == 0:
            skip_lines = 0

        self.df = pd.read_csv(file,skip_blank_lines=True, header=line_num,skiprows=skip_lines,sep="\t")



