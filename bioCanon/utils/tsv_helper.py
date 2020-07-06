import pandas as pd
import os

def parse_tsv_data_to_dict(file,header,logging,header_row=0,index_col=0):
    """

    Parameters
    ----------
    file [str] : Path to TSV file
    header [list] : List of header fields
    logging [logging obj] : Valid logging object
    header_row [int] : row number of TSV file with header
    index_col [int] : column number of TSV file to use as index

    Returns
    -------
    data [dict] : TSV data indexed by index col and named attributes

    """
    if not os.path.isfile(file):
        logging.error("Specified file {} is not found, please check that it exists".format(file))
        return dict()
    if os.path.getsize(file) == 0:
        logging.error("Specified file {} is found but is empty".format(file))
        return dict()

    df = pd.read_csv(file, sep='\t', header=header_row, names=header, index_col=index_col)

    data = {}

    header_len = len(header)

    for index, row in df.iterrows():
        data[index] = {}
        for i in range(0, header_len):
            value = row[i]
            if value == 'nan':
                value = ''
            data[index][header[i]] = row[i]

    return data

def filter_df(df,column_name):
    """
    Accepts a pandas dataframe and a column name and returns a dataframe with all null rows removed from that column
    :param df: [pandas dataframe]
    :param column_name: [str] Name of column in data frame to check for null values
    :return: df [pandas dataframe] where there are no null values in the checked column
    """
    return df[df.column_name.notnull()]