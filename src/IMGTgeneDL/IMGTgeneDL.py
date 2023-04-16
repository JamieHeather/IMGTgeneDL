#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
IMGTupdate.py
A script to download the most up-to-date germline genes from IMGT/GENE-DB
"""


import os
import argparse
import collections as coll
import datetime
import re
import textwrap
import warnings
import sys
from requests import get

# Ensure correct importlib-resources function imported
if sys.version_info < (3, 9):
    import importlib_resources                              # PyPI
else:
    import importlib.resources as importlib_resources       # importlib.resources

__version__ = '0.5.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

warnings.filterwarnings('ignore', message='Unverified HTTPS request')


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="IMGTupdate " + str(__version__) + '\n' +
                    ": Download the most up-to-date TCR germline gene information from IMGT/GENE-DB")

    parser.add_argument('-s', '--in_species', required=False, type=str, default='Homo+sapiens',
                        help='Species to download. Use genus/species, separated with a \'+\'. '
                             'Default = \'Homo+sapiens\'.')
    parser.add_argument('-n', '--common_names', action='store_true', required=False, default=False,
                        help='Flag to specify species by common name.')

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
                        help='Flag to download constant region sequences.')

    parser.add_argument('-L', '--loci', required=False, type=str, default='AB',
                        help='Loci to download, as characters/strings (ABGD). Default = \'AB\'.')

    parser.add_argument('-o', '--out_path', required=False, type=str,
                        help='File name for writing out to, overwriting default auto-generated name.')

    parser.add_argument('-m', '--output_mode', required=False, type=str, default='simple',
                        help='Mode to save data in. Current options are "simple" (default) and "stitchr".')

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
        raise IOError("Failed to download IMGT/GENE-DB raw data. Check connection/input options. \n"
                      "URL = " + url_to_dl)


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
            warnings.warn("Detected " + str(len(invalid)) + " invalid loci characters which will be ignored: " + 
                          ', '.join([x for x in invalid]))

        return valid


def determine_genes(in_args):
    """
    :param in_args: Dict of input CLI arguments
    :return: list of single digit characters denoting gene segments to try to downlooad
    """

    genes_to_dl = []
    if in_args['get_all_regions']:
        genes_to_dl = ['L', 'V', 'D', 'J', 'C']

    else:
        if in_args['get_l']:
            genes_to_dl.append('L')
        if in_args['get_v']:
            genes_to_dl.append('V')
        if in_args['get_d']:
            genes_to_dl.append('D')
        if in_args['get_j']:
            genes_to_dl.append('J')
        if in_args['get_c']:
            genes_to_dl.append('C')

    warnings.warn("Detected " + str(len(genes_to_dl)) + " gene types to download: " + ', '.join(genes_to_dl) + '.')

    return genes_to_dl


def get_specific_items(search_loci, search_genes, species_constants, basic_gene_types, search_species):
    """
    :param search_loci: List of loci (as single digit characters) to search for (e.g. 'ABGD')
    :param search_genes: List of genes (as single digit characters) to search for (e.g. 'LVDJC')
    :param species_constants: Dict of lists of known exon configurations, from get_species_specific_constants
    :param basic_gene_types: Dict of non-constant region codes:IMGT names
    :param search_species: Str of full genus/species to search for (e.g. 'Homo+sapiens')
    :return:
    """

    out_fasta = ''

    # Loop through all relevant input
    for locus in search_loci:

        # Add the species/locus-specific combination of possible known constant configs to the basic gene type list
        gene_types = basic_gene_types
        gene_types['C'] = species_constants['TR' + locus + 'C']

        for gene in search_genes:

            for gene_type in gene_types[gene]:

                # Use these fields to populate the URL...
                # (accounting for the fact that leaders are technically part of the V)
                if gene == 'L':
                    search_gene = 'V'
                else:
                    search_gene = gene

                # Skip D gene search for those loci which don't have any
                if search_gene == 'D' and locus not in ['B', 'D']:
                    continue

                full_gene = 'TR' + locus + search_gene

                # If looking for constant regions, download each exon individually and compile into full sequences
                # (As currently using EX1+EX2...etc attempts fail for constants with exons from different accessions)
                if search_gene == 'C':

                    exons = gene_type.split('+')
                    exon_fasta = ""
                    for exon in exons:
                        dl_url = imgt_url(full_gene, search_species, exon)
                        try:
                            reads = page_scrape(dl_url)
                            if reads:
                                exon_fasta += reads
                            else:
                                warnings.warn("No constant region exons found for: " + '/'.join(
                                    [search_species, full_gene, exon]))

                        except:
                            warnings.warn("Failed to download individual exon IMGT/GENE-DB raw data (for " +
                                          '/'.join([full_gene, search_species, gene_type]) +
                                          ". Check connection/input options")

                    assembled_fasta = concat_constants(exon_fasta, gene_type)
                    out_fasta += assembled_fasta

                # Otherwise try to download the full gene region
                else:
                    dl_url = imgt_url(full_gene, search_species, gene_type)

                    try:
                        # ... and then download that, pulling out the FASTA sequence from amidst the HTML
                        reads = page_scrape(dl_url)
                        if reads:
                            out_fasta += reads

                    except:
                        warnings.warn("Failed to download specific IMGT/GENE-DB raw data for " +
                                      '/'.join([search_species, full_gene, gene_type]) +
                                      ". Check connection/input options. \nURL = " + dl_url)

            # TODO filter out TRAV/DV duplicates (which get downloaded both from TRA and TRD loci calls) at this stage?

    return out_fasta


def imgt_url(four_char_gene, genus_species, imgt_label):
    """
    :param four_char_gene: str of ull gene type descriptor, e.g. TRAC/TRBV/TRGJ etc
    :param genus_species: str of genus/species, separated by a +, e.g. Homo+sapiens
    :param imgt_label: str of IMGT label to download, e.g. V-REGION/EX1+EX2+EX3+EX4 etc
    :return: str of properly concatenated URL to download from
    """
    return 'https://www.imgt.org/genedb/GENElect?query=8.1+' + four_char_gene + \
           '&species=' + genus_species + '&IMGTlabel=' + imgt_label


def page_scrape(download_url):
    """
    :param download_url: URL of IMGT page to try to download FASTA sequences from
    :return: FASTA reads found on the page, if any
    """
    scrape = get(download_url, verify=False).content.decode()

    if '\n>' in scrape:
        fasta_start = scrape.index('\n>')

        if '\n<' in scrape[fasta_start:]:
            fasta_end = scrape[fasta_start:].index('\n<')
            fasta = scrape[fasta_start + 1:fasta_start + fasta_end]
            return fasta

        else:
            return

    else:
        return


def import_common_names(path_to_species):
    """
    :param path_to_species: str path to a file containing a table converting common animal names to genus species
    :return: a dict converting common names to scientific (in IMGT URL format)
    """

    common_names = coll.defaultdict()
    if os.path.isfile(path_to_species):
        with open(path_to_species, 'r') as in_file:
            line_count = 0
            for line in in_file:
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    species_headers = bits
                else:
                    common_names[bits[0]] = bits[1]
                line_count += 1
        return common_names

    else:
        raise IOError("Unable to locate file of species names: " + path_to_species)


def get_species_specific_constants(search_species, default_const_exons, variant_file):
    """
    :param search_species: str of genus/species (separated with a '+' character), e.g. Homo+sapiens
    :param default_const_exons: dict of list of default exon arrangements for each of the loci's constant regions
    :param variant_file: str path to file containing constant regions with non-standard exon arrangements
    :return: A modified version of the gene_type_dict population with the appropriate constant region configurations
    """

    try:
        with open(variant_file, 'r') as in_file:
            line_count = 0

            for line in in_file:
                bits = line.rstrip().split('\t')
                if line_count == 0:
                    headers = bits

                else:
                    # Only look for sequence matches
                    if bits[0] == search_species:
                        default_const_exons[bits[1]].append(bits[2])

                line_count += 1

        return default_const_exons

    except:
        warnings.warn("Failed to read the constant region variant file (" + variant_file + ") in - "
                      "please ensure it's present and correct.")
        return default_const_exons


def output_stitchr_format(cli_args, fasta_str):
    """
    :param cli_args: dict of command line arguments from argparse
    :param fasta_str: str containing all fasta reads concatenated together, as downloaded from IMGT
    :return: the path to the output directory
    """

    # Check output directory
    out_dir = cli_args['common_name'] + '/'  # TODO allow specification of path to output folder?
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # First write out the total FASTA file
    with open(out_dir + 'imgt-data.fasta', 'w') as output_file:
        output_file.write(fasta_str)

    # Then write the individual chain files (note that I include all TRDV genes in the TRA file, as they do recombine)
    out_strs = {'A': "", 'B': "", 'G': "", 'D': ""}
    region_counts = coll.defaultdict(nest_counter)
    j_str = "J gene\tResidue\tConfident?\tMotif\tPosition\n"
    c_str = "C gene\tExons\tStart motif\tStop codon motif\n"

    for record in fasta_str.split('>')[1:]:
        read_header = record[:record.index('\n')]
        read_header_bits = read_header.split('|')
        formatted_seq = record[record.index('\n')+1:]
        read_seq = formatted_seq.replace('\n', '')

        # Note that some IMGT gene sequences contain non-ACGT characters; best just filter the whole reads
        if not dna_check(read_seq):
            warnings.warn("Skipping gene sequence due to non-DNA character detected: " + read_header)
            continue

        # Add an additional field to the FASTA header to provide a consistent gene type label
        read_region = read_header_bits[4]
        if read_region == 'L-PART1+L-PART2':
            new_label = 'LEADER'
        elif read_region == 'V-REGION':
            new_label = 'VARIABLE'
        elif read_region == 'D-REGION':
            new_label = 'DIVERSITY'
        elif read_region == 'J-REGION':
            new_label = 'JOINING'
        elif read_region.startswith('EX'):
            new_label = 'CONSTANT'

            # Account for the fact assembled constant regions seem to acquire an extra space at the end of their headers
            if read_header.endswith(' '):
                read_header = read_header[:-1]

        # Skip partial genes and D-REGIONS, as they're no good for stitching with
        if 'partial' in read_header or 'D-REGION' in read_header:
            continue

        gene_section = read_header_bits[1]
        out_record = '>' + read_header + '~' + new_label + '\n' + formatted_seq

        # First do the straightforward read binning
        # Use nested ifs for the conditionals that need to go in two, but also ensure each only goes in once
        read_locus = gene_section[2]
        out_loci = []
        if read_locus == 'A':
            out_loci.append('A')
            # Joint TRAx/DVy genes can also rearrange with delta
            if '/DV' in gene_section:
                out_loci.append('D')
        elif read_locus == 'B':
            out_loci.append('B')
        elif read_locus == 'G':
            out_loci.append('G')
        elif read_locus == 'D':
            out_loci.append('D')
            # TRDV genes can also rearrange with alpha, so include them there
            if gene_section[3] == 'V':
                out_loci.append('A')

        # Then add that read to the appropriate locus' dict (assuming it hasn't been added already, to prevent dupes)
        for ol in out_loci:
            if out_record not in out_strs[ol]:
                out_strs[ol] += out_record

        # While parsing the reads, we also need to infer J gene conserved CDR3-ending motifs
        if gene_section[3] == 'J':
            j_motifs = determine_j_motifs(read_header_bits, read_seq)
            j_str += '\t'.join(j_motifs) + '\n'

        # Similarly, infer motifs necessary for proper stitchr frame determination
        if gene_section[3] == 'C':
            c_motifs = determine_c_motifs(read_header_bits, read_seq)
            c_str += '\t'.join(c_motifs) + '\n'

        # Count which regions are being found, to ensure we have at least some of each for stitching
        region_counts[read_locus][new_label] += 1

    for loc in out_strs:
        # Only write out loci which were detected...
        if len(out_strs[loc]) > 0:

            # ... and only if they have at least some of each type of read for that locus
            if len(region_counts[loc]) >= 4:
                with open(out_dir + 'TR' + loc + '.fasta', 'w') as output_file:
                    output_file.write(out_strs[loc])
            else:
                warnings.warn("Warning: less than the complete number of sequence types (L/V/J/C) detected for "
                              + loc + " locus; this cannot be stitched, so this locus will be skipped. ")

        else:
            warnings.warn("Warning: no sequences detected for " + loc + " locus. ")

    # Write out J/C region motif data
    with open(out_dir + 'J-region-motifs.tsv', 'w') as output_file:
        output_file.write(j_str)

    with open(out_dir + 'C-region-motifs.tsv', 'w') as output_file:
        output_file.write(c_str)

    # Finally generate the brief summary log file
    with open(out_dir + 'data-production-date.tsv', 'w') as output_file:
        out_str = '\n'.join([
            'imgt-data.fasta_last_modified\t' + today(),
            'script_used\tIMGTgeneDL',
            'last_run\t' + today(),
            'version_used\t' + __version__,
            'imgt_genedb_release\t' + get_release(base_url)
        ])
        output_file.write(out_str)

    return out_dir


def determine_c_motifs(gene_header_bits, gene_seq):
    """
    We need recognisable in-frame amino acid sequences at the start of each constant region, and after any stop codons,
    to allow stitchr to correctly infer frame and trim coding sequences appropriately.
    :param gene_header_bits: list of IMGT FASTA header sections, having split on pipe ('|') character
    :param gene_seq: str of FASTA sequence of gene
    :return: List of [gene, exon config, starting motif, stop codon motif (if any)]
    """

    func, seq_type = gene_header_bits[3:5]

    # Check to see whether the IMGT C gene is cDNA only (as those include an additional 5' nt from the J)
    if func.startswith('(') and func.endswith(')'):
        translation = translate(gene_seq[3:])
    else:
        translation = translate(gene_seq[2:])

    c = gene_header_bits[1]
    start_motif = translation[:10]
    stop_motif = ''
    if '*' in translation:
        stop_site = translation.index('*')
        stop_motif = translation[stop_site:stop_site + 11]

    return [c, seq_type, start_motif, stop_motif]


def determine_j_motifs(gene_header_bits, gene_seq):
    """
    Aims to automatically infer what the conserved CDR3-ending motif is for any given J gene
    Note that the core of this code is taken from autoDCR/generate-tag-files.py
    See repo here: https://github.com/JamieHeather/autoDCR
    :param gene_header_bits: list of IMGT FASTA header sections, having split on pipe ('|') character
    :param gene_seq: str of FASTA sequence of gene
    :return: list of [gene, first residue (F), confidence (Y/N/0), 4-digit motif (FGFG), position (-11)]
    """

    poss_motifs = ['FG.G', 'WG.G', 'CG.G', 'F..G', 'FG..', 'LG.G']
    j = gene_header_bits[1]
    func = gene_header_bits[3]

    # Define the expected motif start site position(s)
    expected_starts = {
        'TRAJ': [-34],
        'TRBJ': [-31],
        'TRGJ': [-37, -31],
        'TRDJ': [-34]
    }

    found = False
    for start in expected_starts[j[:4]]:
        # Check to see whether the IMGT J gene is cDNA only (as those omit the last nucleotide)
        if func.startswith('(') and func.endswith(')'):
            start += 1
            modifier = 1
        else:
            modifier = 0

        base = translate(gene_seq[start:start + 12])

        # Get AA seq (NB: end, not start of the J needs to be in frame, after deleting the base that splices to C)
        modulo = (len(gene_seq) - 1 + modifier) % 3
        translation = translate(gene_seq[modulo:])

        # Record if we get a hit at the appropriate side
        if [y[0] for y in [re.findall(x, base) for x in poss_motifs] if y]:
            found = True
            motif = base
            confidence = 'Y'
            break

    # If nothing found after looking at conserved sites (in full-length sequences), then look across whole sequence
    if not found:
        search = list(set([y[0] for y in [re.findall(x, translation) for x in poss_motifs] if y]))
        if len(search) == 1:
            seq_search = (translation.index(search[0]) * 3) + modulo
            start = seq_search - len(gene_seq)
            motif = search[0]
            confidence = 'N'

        #  Otherwise there's just no motifs to be found!
        else:
            motif = base
            confidence = '0'

    position = str(int((start + 1 - modifier) / 3))
    return [j, motif[0], confidence, motif, position]


def dna_check(possible_dna):
    """
    :param possible_dna: A sequence that may or may not be a plausible DNA (translatable!) sequence
    :return: True/False
    """

    return set(possible_dna.upper()).issubset({'A', 'C', 'G', 'T'})


def nest():
    """
    :return: an instance of a list defaultdict (for the creation of nested ddicts)
    """
    return coll.defaultdict(list)


def nest_counter():
    """
    Create nested counters
    """
    return coll.Counter()


def translate(seq):
    """
    :param seq: Nucleotide sequence
    :return: Translated nucleotide sequence, with overhanging 3' (J) residues trimmed
    """
    protein = ""
    # Trim sequence length to a multiple of 3 (which in this situation should just be removing the terminal J residue)
    difference = len(seq) % 3
    if difference != 0:
        seq = seq[:-difference]
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].upper()
        protein += codons[codon]
    return protein


def concat_constants(exon_reads, exon_arrangement):
    """
    :param exon_reads: str of all extracted reads for the individual exons
    :param exon_arrangement: str of expected configuration of the exons in the final spliced constant region
    :return: str of fasta of the assembled full spliced constant regions
    """

    # First read in all the downloaded individual exons to gene-specific dicts
    gene_bits = coll.defaultdict(nest)
    header_bits = coll.defaultdict(nest)
    reads = exon_reads.split('>')
    for r in reads:
        if r:
            header_break = r.index('\n')
            header = r[:header_break]
            seq = r[header_break + 1:].replace('\n', '').replace('\r', '')
            bits = header.split('|')
            gn = bits[1]
            ex = bits[4]
            if ex not in gene_bits[gn]:
                gene_bits[gn][ex] = seq
                header_bits[gn][ex] = bits
            else:
                raise IOError("Exon clash detected in " + gn + "/" + ex)

    # Then go through each gene and try to assemble the desired configuration of exons
    out_fasta = ""
    goal = exon_arrangement.split('+')
    for g in gene_bits:
        incomplete = False
        gene_seq = ""
        for x in goal:
            if x in gene_bits[g]:
                gene_seq += gene_bits[g][x]
            else:
                # print("Potential missing exon: " + '/'.join([g, x]))
                incomplete = True

        if not incomplete:
            # If it made it this far, then the complete sequence has been made - so we can make the appropriate header
            accessions = ','.join(list(set([header_bits[g][y][0] for y in header_bits[g]])))
            functionality = ','.join(list(set([header_bits[g][y][3] for y in header_bits[g]])))
            rebuilt_exons = '+'.join([header_bits[g][y][4] for y in header_bits[g]])

            if rebuilt_exons != exon_arrangement:
                raise IOError("Exon configuration not correct for " + g + ": " + rebuilt_exons)
            len_str = str(len(gene_seq))

            # We also need to tweak the gene identifier for non-standard exon configurations, to prevent duplicates
            suffix = ''
            if exon_arrangement not in [x[0] for x in constant_regions_exons.values()]:
                suffix = exon_arrangement
                for substr in ['EX1', 'EX2', 'EX3', 'EX4', '+', 'UTR']:
                    suffix = suffix.replace(substr, '')
                suffix = '_' + suffix

            out_header = '|'.join([accessions, g + suffix, header_bits[g]['EX1'][2], functionality, exon_arrangement,
                                   '?', len_str + ' nt', '?', ' ', ' ', ' ', ' ', len_str + '+0=' + len_str, ' ', ' '])
            out_fasta += fastafy(out_header, gene_seq)

    return out_fasta


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

codons = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
          'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
          'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
          'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
          'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
          'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
          'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
          'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
          'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
          'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
          'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
          'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
          'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
          'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
          'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
          'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}


# Note that this dict provides a list of default regions to download
# Constant regions are not featured here, but calculated on the fly due to their changing between loci/species
default_gene_types = {
    'L': ['L-PART1+L-PART2'],
    'V': ['V-REGION'],
    'D': ['D-REGION'],
    'J': ['J-REGION']
}

# Specifies the expected (minimal) exons required to produce a coding constant region
# Ignores less conserved/duplicated TRGC exon configurations, set via the c-region-variant-configs.tsv file
constant_regions_exons = {
    'TRAC': ['EX1+EX2+EX3+EX4UTR'],
    'TRBC': ['EX1+EX2+EX3+EX4'],
    'TRGC': ['EX1+EX2+EX3'],
    'TRDC': ['EX1+EX2+EX3+EX4UTR']
}


def process_input_args(cli_args, default_species):
    """
    :param cli_args: dict of command line arguments from argparse
    :return: modified dictionary having sorted various provisions
    """
    # First sort species information

    # If the common species name has been provided, check that it is present in the default list ...
    species_conversion = import_common_names(default_species)

    named = False
    # If common name used
    if 'common_names' in cli_args:  # Belt and braces approach to testing if variable present, works on more platforms
        if cli_args['common_names']:

            if cli_args['in_species'].upper() in species_conversion:
                cli_args['species'] = species_conversion[cli_args['in_species'].upper()]
                cli_args['common_name'] = cli_args['in_species'].upper()
            else:
                raise IOError("Common name input mode selected (-n) but species '" + cli_args['species'].upper() +
                              "' not in species file. Please check naming/file (" + default_species + ").")

            named = True

    # Otherwise, scientific name must be used
    if not named:
        if '+' not in cli_args['in_species']:
            raise IOError("Error in input species name: '+' characters expected (for non-common names).")
        else:
            if cli_args['in_species'].upper() in [x.upper() for x in species_conversion.values()]:
                search = [x for x in species_conversion if
                          species_conversion[x].upper() == cli_args['in_species'].upper()][0]
                cli_args['species'] = species_conversion[search]
                cli_args['common_name'] = search

    # Determine what users want to download
    # if cli_args['output_mode'] == 'simple':

    # Stitchr mode requires all loci/all regions be downloaded
    if cli_args['output_mode'] == 'stitchr':
        print("Processing IMGT data into Stitchr data format.")
        cli_args['get_all_regions'] = True
        cli_args['loci'] = 'TR'

    return cli_args


def main():

    # Sort input/output file details
    pkg_files = importlib_resources.files("IMGTgeneDL")
    c_region_variants_file = str(pkg_files / 'c-region-variant-configs.tsv')
    default_species_path = str(pkg_files / 'species.tsv')
    print(str(pkg_files))
    print(importlib_resources.files("IMGTgeneDL"))
    input_args = process_input_args(vars(args()), default_species_path)

    with warnings.catch_warnings(record=True) as warnings_list:
        warnings.simplefilter("always")

        if input_args['get_all']:
            url, out_suffix = get_url(input_args, base_url)
            out_path = name_out_file(input_args, out_suffix)
            download_genedb(url, out_path)

        else:

            # Determine which loci and gene segments to get
            constant_regions = get_species_specific_constants(input_args['species'],
                                                              constant_regions_exons,
                                                              c_region_variants_file)

            loci = determine_loci(input_args)
            genes = determine_genes(input_args)
            out_path = name_out_file(input_args, [input_args['species'], ''.join(loci), ''.join(genes)])

            # Then scrape IMGT for those sections
            compiled_sequences = get_specific_items(loci, genes, constant_regions,
                                                    default_gene_types, input_args['species'])

            # Determine the appropriate output format, with the default being a 'simple' output list
            if input_args['output_mode'] == 'simple':
                with open(out_path, 'w') as out_file:
                    out_file.write(compiled_sequences)

            # Alternatively output in the appropriate format for Stitchr
            elif input_args['output_mode'] == 'stitchr':
                stitchr_dir = output_stitchr_format(input_args, compiled_sequences)

            else:
                raise IOError("Unknown output mode specified: " + input_args['output_mode'])

    # Finally output a warnings log
    warnings_txt = '\n'.join([str(warnings_list[x].message) for x in range(len(warnings_list)) if
                              'HTTPS' not in str(warnings_list[x].message) and
                              'certificate' not in str(warnings_list[x].message)])
    warnings_file = 'IMGTgeneDLwarnings.txt'
    if input_args['output_mode'] == 'stitchr':
        warnings_file = stitchr_dir + warnings_file

    with open(warnings_file, 'w') as out_file:
        out_file.write(warnings_txt)
