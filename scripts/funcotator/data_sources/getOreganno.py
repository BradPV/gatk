#!/usr/bin/env python

########################################################################
# Imports:

import csv
import sys
import urllib2
from GenericTsvReader import GenericTsvReader

########################################################################
# Constants:

FILE_URL = 'http://www.oreganno.org/dump/ORegAnno_Combined_2016.01.19.tsv'
RAW_FILE_NAME = FILE_URL.split('/')[-1]

OUTPUT_HEADERS = ['Build', 'Chr', 'Start', 'End', 'ID', 'Values']

OUT_HG19_FILE_NAME_TEMPLATE = 'oreganno_%s.HG19.tsv'
OUT_HG38_FILE_NAME_TEMPLATE = 'oreganno_%s.HG38.tsv'

OREGANNO_KEYS = ['Outcome', 'Type', 'Gene_Symbol', 'Gene_ID', 'Gene_Source', 'Regulatory_Element_Symbol',
                'Regulatory_Element_ID', 'Regulatory_Element_Source', 'dbSNP_ID', 'PMID', 'Dataset']

VALUES_DELIMITER = '|'
NON_VALUE = "N/A"

########################################################################
# Functions:

FLUSH = sys.stdout.flush


def get_values_data_from_row_dict(r):
    values = []

    for i in xrange(0, len(OREGANNO_KEYS)):
        key = OREGANNO_KEYS[i]
        if r[key] != NON_VALUE:
            values.append(key + "=" + r[key])

    return VALUES_DELIMITER.join(values)


def get_file_from_url(url, file_name):
    """Taken from: https://stackoverflow.com/a/22776"""
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    print "Downloading: %s Bytes: %s" % (file_name, file_size)

    file_size_dl = 0
    block_sz = 8192
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8) * (len(status) + 1)
        print status,

    f.close()

########################################################################
# Main:

if __name__ == '__main__':

    print 'Downloading Oreganno file:', RAW_FILE_NAME, '...',
    FLUSH()
    # Download the Raw file:
    get_file_from_url(FILE_URL, RAW_FILE_NAME)
    print 'DONE!'
    FLUSH()

    # Get the version of the file:
    file_version = RAW_FILE_NAME.replace('ORegAnno_Combined_', '').replace('.tsv', '').replace('.', '')

    print 'Creating output files ...',
    FLUSH()
    # Set up the output files:
    out_hg19_tsv_writer = csv.DictWriter(file(OUT_HG19_FILE_NAME_TEMPLATE % file_version, 'w'), OUTPUT_HEADERS,
                                         delimiter='\t', lineterminator="\n")
    out_hg19_tsv_writer.fieldnames = OUTPUT_HEADERS
    out_hg19_tsv_writer.writeheader()
    out_hg38_tsv_writer = csv.DictWriter(file(OUT_HG38_FILE_NAME_TEMPLATE % file_version, 'w'), OUTPUT_HEADERS,
                                         delimiter='\t', lineterminator="\n")
    out_hg38_tsv_writer.fieldnames = OUTPUT_HEADERS
    out_hg38_tsv_writer.writeheader()
    print 'DONE!'
    FLUSH()

    # Now that we have the file, go through it and reformat it:
    tsvReader = GenericTsvReader(RAW_FILE_NAME)
    headers = tsvReader.getFieldNames()
    print('Found headers (input): ' + str(headers))
    FLUSH()

    print 'Iterating through input file ...',
    FLUSH()
    # Go through our file:
    for lineNum, line in enumerate(tsvReader):

        if line['Species'].lower() != 'homo sapiens':
            print '    Ignoring non-human record on line:', lineNum, " : ", line['Species'].lower()
            FLUSH()
            continue

        # Get the trivial fields here:
        row = dict()
        row['Build'] = line['Build']
        row['Chr'] = line['Chr']
        row['Start'] = line['Start']
        row['End'] = line['End']
        row['ID'] = line['ORegAnno_ID']

        # Get the values field from our helper method:
        row['Values'] = get_values_data_from_row_dict(line)

        if line['Build'].lower() == 'hg19':
            out_hg19_tsv_writer.writerow(row)
        elif line['Build'].lower() == 'hg38':
            out_hg38_tsv_writer.writerow(row)
        else:
            print '    Ignoring incompatible build record on line:', lineNum, " : ", line['Build'].lower()
            FLUSH()
    print 'DONE!'
    FLUSH()


