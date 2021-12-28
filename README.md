# IMGTgeneDL

## 0.1.0
##### Jamie Heather | CCR @ MGH | 2021

This script provides an alternative way to access TCR and IG genes stored in [IMGT/GENE-DB](http://www.imgt.org/genedb/). It's primarily designed for downloading human/mouse TCRs, but it's readily adaptable to other species/loci.

## Usage

This script is tested on python >= 3.6, not requiring any non-standard packages.

### Download specific gene types

The primary way this script is intended to be used is to tell it the species, loci, and sequence types that you want to download. This will be downloaded and saved to a file named with the date the [IMGT release](http://www.imgt.org/download/GENE-DB/RELEASE) used, and details of the combination of parameters searched for (unless overriden with the `-o / --out_path` flag).

#### Species

The script must be run on single species at a time, given via the `-s / --species` flag as a full genus species with a '+' symbol in place of the space. E.g.:

* `-s Homo+sapiens`
* `-s Mus+musculus`

Note that it doesn't seem that the IMGT URL interface will accept either genus or species alone, and it can be particular about formatting so maintaining proper case is advised. However it seems to download sub-species (e.g.\ searching for Mus+musculus will return all covered strains).

#### Loci

This script is currently configured to download the four common TCR loci:

* A / TRA / alpha
* B / TRB / beta
* G / TRG / gamma
* D / TRD / delta

These must be provided to the script using the `-L / --loci` flag, giving it the desired loci as a string of characters, e.g. `-L AB` to just download alpha and beta sequences, or `-L G` to just download gamma. Alternatively `-L TR` will simply download all four chains' sequences (equivalent to `-L ABGD`).

#### Sequence types

This script is designed to help in the aid in the analysis of typical expressed repertoires, and thus is configured by default to download the relevant parts of the loci that end up involved in expressed transcripts. The download of each is achieved using specific flags:

* `-l / --get_l`: download leader sequences
* `-v / --get_v`: download V sequences
* `-d / --get_d`: download D sequences
* `-j / --get_j`: download J sequences
* `-c / --get_c`: download constant region sequences

Note that these can be combined, e.g. `-vdj` will just download the V, D, and J gene sequences. Alternatively users can apply the `-r / --get_all_regions` flag to just download all of these regions (equivalent to `-lvdjc`).

#### Examples

The following is the basic command to download all relevant human sequences for all chains:
```
python3 IMGTgeneDL.py -s Homo+sapiens -L TR -r
```

While this is a command to just download delta chain J genes from mice:
```
python3 -i IMGTgeneDL.py -s Mus+musculus -j -L D
```

### Download whole database

If no locus and gene type flags are used, or if the `-a / --get_all` flag is used, then the script will just download [the whole of GENE-DB](http://www.imgt.org/download/GENE-DB/) - all species, all genes, all loci. By default this downloads the [ungapped nucleotide file, with all pseudogenes](http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP), and saves this to a file named with the date and the [IMGT release](http://www.imgt.org/download/GENE-DB/RELEASE) used. This can be changed using the following flags:

* `-gap / --gapped`: downloads the gapped FASTA instead of the ungapped
* `-ifp / --in_frame_p`: downloads the 'inframeP' FASTA instead of the 'allP'
* `-o / --out_path`: as above, sets the path to a specific file if you don't wish to use the automatic names or save in the same directory

## Notes

### Downloaded regions 

The architecture of the TCR loci differs a little between the genes and across species, which the IMGT nomenclature has specific terms to cope with. However the URL based searching this script does requires provision exact **exon names** for each species. This script assumes generic defaults, but these can be overriden by providing specific details in the tab-delimited `region-overrides.tsv` file. This allows users to override which fields are downloaded, or download additional fields, by adding an entry for the relevant gene/species combination and filling in the final 'Field(s)' comma-delimited field with the IMGT labels to be downloaded.

The most relevant place this comes in to place is in the constant regions, which have differing numbers and names of exons. The relevant differences for these loci/species are that exon 4 of the alpha and delta chains is an UTR, while gamma chains lack a fourth exon and have duplicated exon 2 variants. If users wish to run the script to download specific sequences including constant regions for species *other than humans or mice* they will need to edit this document appropriately first.

The other default IMGT labels downloaded are:
* L-PART1+L-PART2 for leader sequences
* V-/D-/J-REGION for V/D/J genes

### IMGT FASTA headers

The IMGT header FASTA fields (as reported in the output of [GENE-DB](http://www.imgt.org/genedb/)) are:

```
The FASTA header contains 15 fields separated by '|':

1. IMGT/LIGM-DB accession number(s)
2. IMGT gene and allele name
3. species
4. IMGT allele functionality
5. exon(s), region name(s), or extracted label(s)
6. start and end positions in the IMGT/LIGM-DB accession number(s)
7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
8. codon start, or 'NR' (not relevant) for non coding labels
9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
12. number of amino acids (AA): this field indicates that the sequence is in amino acids
13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
14. partial (if it is)
15. reverse complementary (if it is)
```

##### Disclaimer

I am not affiliated with IMGT, and this tool is only shared as a way to increase the utility of their platform. Please TCR responsibly.
