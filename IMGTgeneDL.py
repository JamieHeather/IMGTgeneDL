#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
IMGTupdate.py
A script to download the most up-to-date germline genes from IMGT/GENE-DB
"""


import os
import argparse
import datetime
import textwrap
from requests import get

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="IMGTupdate " + str(__version__) + '\n' +
                    ": Download the most up-to-date TCR germline gene information from IMGT/GENE-DB")

    parser.add_argument('-s', '--species', required=False, type=str, default='Homo+sapiens',
                        help='Species to download. Use genus/species, separated with a \'+\'. '
                             'Default = \'Homo+sapiens\'.')

    # Fields related to downloading specific fields
    parser.add_argument('-r', '--get_all_regions', action='store_true', required=False, default=False,
                        help='Flag to download all regions (leader, V, D, J, constant) for the specific loci.'
                             '\nNB: overrides all -l/-v/-d/-j/-c options.')
    parser.add_argument('-l', '--get_l', action='store_true', required=False, default=False,
                        help='Flag to download leader (L-PART1+L-PART2) sequences.')
    parser.add_argument('-v', '--get_v', action='store_true', required=False, default=False,
                        help='Flag to download variable gene (V-REGION) sequences.')
    parser.add_argument('-d', '--get_d', action='store_true', required=False, default=False,
                        help='Flag to download diversity gene (D-REGION) sequences.')
    parser.add_argument('-j', '--get_j', action='store_true', required=False, default=False,
                        help='Flag to download joining gene (J-REGION) sequences.')
    parser.add_argument('-c', '--get_c', action='store_true', required=False, default=False,
                        help='Flag to download constant region (EX1/EX2/EX3/EX4) sequences.')

    parser.add_argument('-L', '--loci', required=False, type=str, default='AB',
                        help='Loci to download, as characters/strings (ABGD). Default = \'AB\'.')

    parser.add_argument('-o', '--out_path', required=False, type=str,
                        help='File name for writing out to, overwriting default auto-generated name.')

    # Fields related to downloading the whole database
    parser.add_argument('-a', '--get_all', action='store_true', required=False, default=False,
                        help='Disregard all specific download options and just grab whole GENE-DB database.'
                             '\nNB: Only contains V/D/J/C sequence (i.e. no leader sequences).')

    parser.add_argument('-gap', '--gapped', action='store_true', required=False, default=False,
                        help='Flag to download gapped genes if using the -a flag, rather than ungapped (default).')

    parser.add_argument('-ifp', '--in_frame_p', action='store_true', required=False, default=False,
                        help='Optional flag to download in-frame pseudogenes if using the -a flag,'
                             ' rather than all (default).')

    return parser.parse_args()


def download_genedb(url_to_dl, fasta_file):
    """
    :param url_to_dl: URL of requested IMGT-GENE/DB resource
    :param fasta_file: path to fasta file to save locally
    :return: nothing, just save the contents of the URL locally
    """

    try:
        with open(fasta_file, "wb") as file:
            response = get(url_to_dl, verify=False)
            file.write(response.content)

    except:
        raise IOError("Failed to download IMGT/GENE-DB raw data. Check connection/input options")


def readfa(fp):
    """
    readfq(file):Heng Li's Python implementation of his readfq function
    https://github.com/lh3/readfq/blob/master/readfq.py
    :param fp: opened file containing fastq or fasta reads
    :yield: read id, read sequence, and (where available) read quality scores
    """

    last = None  # this is a buffer keeping the last unprocessed line
    while True:

        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] == '>':  # fasta header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:], [], None  # This version takes the whole line (post '>')
        for l in fp:  # read the sequence
            if l[0] == '>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs)  # yield a fasta record
            if not last:
                break


def today():
    """
    :return: Today's date, in ISO format
    """
    return datetime.datetime.today().date().isoformat()


def name_out_file(in_args, other_fields):
    """
    :param in_args: tidied dictionary of arguments from the argparse CLI, for use in 'all' mode
    :param other_fields: list of other strings to be included in filename for 'all/specific' modes
    :return: str of name (and potentially path) of outfile
    """

    # See if there's a user-provided out file name
    if in_args['out_path']:
        out_file_name = in_args['out_path']

        # If a path is provided, check it exists/is valid
        if '/' in out_file_name:
            if out_file_name.endswith('/'):
                raise IOError("Provided path does not include a file name.")
            if not os.path.isdir(out_file_name[:out_file_name.rfind('/')]):
                raise IOError("Provided path is not found.")

        # Then check the provided name has the right extension
        if not out_file_name.endswith('.fasta') and not out_file_name.endswith('.fa'):
            out_file_name += '.fasta'

    # Otherwise infer a descriptive name from the parameters of the data extracted
    else:
        release = get_release(base_url)

        out_file_name = '_'.join([today(), release] + other_fields) + '.fasta'

    return out_file_name


def get_url(in_args, url_stem):
    """
    :param in_args: CLI arguments
    :param url_stem: base URL for all the IMGT GENE-DB resources
    :return: 2 str - the appropriate IMGT FASTA URL & the suffix denoting what's in the downloaded file
    """

    if not in_args['gapped'] and not in_args['in_frame_p']:
        return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP", \
                ['IMGT-GENEDB', 'ungapped', 'allP']

    elif not in_args['gapped'] and in_args['in_frame_p']:
        return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP", \
                ['IMGT-GENEDB', 'ungapped', 'inframeP']

    elif in_args['gapped'] and in_args['in_frame_p']:
        return url_stem + "IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP", \
                ['IMGT-GENEDB', 'gapped', 'inframeP']

    elif in_args['gapped'] and not in_args['in_frame_p']:
        raise IOError("IMGT/GENE-DB does not offer a gapped nucleotide FASTA with in frame pseudogenes.")


def get_release(url_stem):
    """
    :param url_stem:
    :return: the release number taken from the IGMT website
    """

    try:
        return get(url_stem + 'RELEASE', verify=False).content.decode()

    except:
        raise IOError("Unable to download GENE-DB release number. Please check connection.")


def fastafy(fasta_header, fasta_seq):
    """
    :param fasta_header: Gene symbol, extracted from the read id
    :param fasta_seq: Total protein primary sequence, extracted from input FASTA/generated by in silico splicing
    :return: An output-compatible FASTA entry ready for writing to file
    """
    return ">" + fasta_header + "\n" + textwrap.fill(fasta_seq, 60) + "\n"


def determine_loci(in_args):
    """
    :param in_args: Dict of input CLI arguments
    :return: list of single digit characters denoting (valid) loci to search for
    """

    if in_args['loci'] == 'TR':
        print("All TR loci option detected.")
        return list(loci_ids.keys())

    else:
        input_loci = [x.upper() for x in in_args['loci']]
        input_loci.sort()
        valid = [x for x in loci_ids if x in input_loci]
        invalid = [x for x in input_loci if x not in loci_ids]
        print("Detected " + str(len(valid)) + " valid loci to download: " + ', '.join([loci_ids[x] for x in valid]))
        if invalid:
            print("Detected " + str(len(invalid)) + " invalid loci characters which will be ignored: " +
                  ', '.join([x for x in invalid]))

        return valid
    # TODO add this info to a log?


def determine_genes(in_args):
    """
    :param in_args: Dict of input CLI arguments
    :return: list of single digit characters denoting gene segments to try to downlooad
    """

    genes_to_dl = []
    if in_args['get_all_regions']:
        genes_to_dl = ['L', 'V', 'D', 'J', 'C']

    else:
        if input_args['get_l']:
            genes_to_dl.append('L')
        if input_args['get_v']:
            genes_to_dl.append('V')
        if input_args['get_d']:
            genes_to_dl.append('D')
        if input_args['get_j']:
            genes_to_dl.append('J')
        if input_args['get_c']:
            genes_to_dl.append('C')

    print("Detected " + str(len(genes_to_dl)) + " gene types to download: "
          + ', '.join(['|'.join(gene_types[x]) for x in genes_to_dl]))
    # TODO add this info to a log?

    return genes_to_dl


def get_specific_items(search_loci, search_genes, search_species):
    """
    :param search_loci: List of loci (as single digit characters) to search for (e.g. 'ABGD')
    :param search_genes: List of genes (as single digit characters) to search for (e.g. 'LVDJC')
    :param search_species: Str of full genus/species to search for (e.g. 'Homo+sapiens')
    :return:
    """

    out_fasta = ''

    # Loop through all relevant input
    for locus in search_loci:
        specific_gene_types = get_specific_gene_types(search_species, locus, gene_types)

        for gene in search_genes:
            for gene_type in specific_gene_types[gene]:

                # Use these fields to populate the URL...
                # (accounting for the fact that leaders are technically part of the V)
                if gene == 'L':
                    search_gene = 'V'
                else:
                    search_gene = gene

                # Skip D gene search for those loci which don't have any
                if search_gene == 'D' and locus not in ['B', 'D']:
                    continue

                dl_url = 'http://www.imgt.org/genedb/GENElect?query=8.1+' + 'TR' + locus + search_gene + \
                         '&species=' + search_species + \
                         '&IMGTlabel=' + gene_type

                try:
                    # ... and then download that, pulling out the FASTA sequence from amidst the HTML
                    page_scrape = get(dl_url, verify=False).content.decode()
                    fasta_start = page_scrape.index('\n>')
                    fasta_end = page_scrape[fasta_start:].index('\n<')
                    fasta = page_scrape[fasta_start + 1:fasta_start + fasta_end]
                    out_fasta += fasta

                except:
                    raise IOError("Failed to download specific IMGT/GENE-DB raw data. Check connection/input options")

    return out_fasta


def get_specific_gene_types(search_species, search_locus, gene_type_dict):
    """
    :param search_species: Str of genus/species (separated with a '+' character), e.g. Homo+sapiens
    :param search_locus: Str of specific locus to search for, e.g. one of 'A/B/G/D'
    :param gene_type_dict: Dict of the default gene type regions
    :return: A modified version of the gene_type_dict
    """

    try:
        with open('region-overrides.tsv', 'r') as in_file:
            line_count = 0

            for line in in_file:
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    headers = bits
                else:
                    # Only look for sequence matches
                    if bits[0] == search_species:
                        if bits[1][2] == search_locus:
                            gene_type_dict[bits[1][3]] = bits[2].split(',')

                line_count += 1
        return gene_type_dict

    except:
        raise IOError("Failed to read the region-overrides.tsv file in - "
                      "please ensure it's present and correct.")
        # TODO allow users to specific directory of override file?


base_url = "https://www.imgt.org/download/GENE-DB/"


# The contents of the IMGT fasta headers that users might wish to filter on
imgt_fields = {
    'accession': 0,
    'gene': 1,
    'species': 2,
    'functionality': 3,
    'length': 6
}

loci_ids = {
    'A': 'alpha',
    'B': 'beta',
    'G': 'gamma',
    'D': 'delta'
}

# Note that this dict provides a list of *default* regions to download
# These can be overriden by changing the 'region-overrides.tsv' file
gene_types = {
    'L': ['L-PART1+L-PART2'],
    'V': ['V-REGION'],
    'D': ['D-REGION'],
    'J': ['J-REGION'],
    'C': ['EX1', 'EX2', 'EX3', 'EX4']
}


if __name__ == '__main__':

    # Sort input/output file details
    input_args = vars(args())

    # Determine what users want to download
    # If no specific sequence types are being requested, or if -a flag set, then just download the whole file
    if input_args['get_all'] or \
            (not input_args['get_all_regions'] and
             not input_args['get_l'] and
             not input_args['get_v'] and
             not input_args['get_d'] and
             not input_args['get_j'] and
             not input_args['get_c']):

        url, out_suffix = get_url(input_args, base_url)
        out_path = name_out_file(input_args, out_suffix)
        download_genedb(url, out_path)

    else:

        # Determine which loci and gene segments to get
        loci = determine_loci(input_args)
        genes = determine_genes(input_args)
        species = input_args['species'].replace(' ', '+')
        out_path = name_out_file(input_args, [species, ''.join(loci), ''.join(genes)])

        # Then scrape IMGT for those sections
        compiled_sequences = get_specific_items(loci, genes, input_args['species'])

        with open(out_path, 'w') as out_file:
            out_file.write(compiled_sequences)

    # TODO add log?
    # TODO remove verify=False when IMGT sorts out their SSL
