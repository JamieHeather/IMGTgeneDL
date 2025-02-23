# IMGTgeneDL

## 0.6.1
##### Jamie Heather | KF-CCR @ MGH | Last updated: 2025-02

This tool provides an alternative way to access TCR and IG genes stored in [IMGT/GENE-DB](http://www.imgt.org/genedb/) via the command-line, aiming to:

* Make it easier to incorporate IMGT sequences into analysis pipelines
* Keep those sequences updated
* Maintain information as to the provenance of the data

It was mostly built to download human and mouse TCR sequences, with so that's what is most tested, but it's readily adaptable to other species and loci.

## Usage

This version has been tested on Python 3.11 and 3.12. It requires Python >3.9.

### Installation

```bash
pip install IMGTgeneDL
```

### Download specific gene types

The primary way this script is intended to be used is to tell it the species, loci, and sequence types that you want to download. This will be downloaded and saved to a file named with the date the [IMGT release](http://www.imgt.org/download/GENE-DB/RELEASE) used, and details of the combination of parameters searched for (unless overriden with the `-o / --out_path` flag).

#### Species

The script must be run on single species at a time, given via the `-s / --species` flag as a full genus species with a '+' symbol in place of the space. E.g.:

* `-s Homo+sapiens`
* `-s Mus+musculus`

Note that it doesn't seem that the IMGT URL interface will accept either genus or species alone, and it can be particular about formatting so maintaining proper case is advised. However it seems to download sub-species (e.g.\ searching for Mus+musculus will return all covered strains).

Alternatively, users can use the `-n / --common_names` flag, and provide a common name. This only works for species that are included in the provided `species.tsv` file, which maps a common name into a scientific one for all the species for which sufficient TCR gene data exists in IMGT at the time of writing. This can be specified like so:

`... -n -s human ...`
`... -n -s SHEEP ...`

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

##### Examples

The following is the basic command to download all relevant human sequences for all chains:
```bash
IMGTgeneDL -s Homo+sapiens -L TR -r
```

While this is a command to just download delta chain J genes from mice:
```bash
IMGTgeneDL -n -s mouse -j -L D
```

#### Downloading any IMGT sequence type

IMGT genes are annotated with a large number of other sequence features that users may wish to download, in addition to the basic coding features the flags described above. To download these users can make use of the `-d / --get_d` command line input, with multiple regions able to be requested, separated by commas, with any apostrophes escaped with a backslash. As all fields need to relate to a V/D/J/C gene, these must also be specified, after an '=' symbol. For example:

```bash
# Download TRBJ recombination signals
IMGTgeneDL -L B -x J-RS=J

# Download TRA 'V-GENE' and leader part 1 regions
IMGTgeneDL -L A -x V-GENE=V,L-PART1=V

# Download TRD 5' UTRs
IMGTgeneDL -L D -x 5\'UTR=V
```

Note that output files produced in this manner will not be labelled with every user-specified gene type - instead output files will be labelled 'X1', or 'X1X2...' and so on, depending on how many fields were requested. 

In order to allow for IMGT introducing or changing field names none are hard-coded into `IMGTgeneDL`. Field names can be identified by looking at the IMGT documentation, or by running [a query via the IMGT/GENE-DB web interface](https://www.imgt.org/genedb/). Please ensure you're getting the field name correct, as `IMGTgeneDL` won't give an error message if an incorrect name is returned (as it cannot distinguish that from a species/locus combination that just lacks data on IMGT). 

### Download whole database

If no locus and gene type flags are used, or if the `-a / --get_all` flag is used, then the script will just download [the whole of GENE-DB](http://www.imgt.org/download/GENE-DB/) - all species, all genes, all loci. By default this downloads the [ungapped nucleotide file, with all pseudogenes](http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP), and saves this to a file named with the date and the [IMGT release](http://www.imgt.org/download/GENE-DB/RELEASE) used. This can be changed using the following flags:

* `-gap / --gapped`: downloads the gapped FASTA instead of the ungapped
* `-ifp / --in_frame_p`: downloads the 'inframeP' FASTA instead of the 'allP'
* `-o / --out_path`: as above, sets the path to a specific file if you don't wish to use the automatic names or save in the same directory

### Output modes

This script has the option to specify different types of mode, selected via the `-m / --output_mode` flag. Currently there are two options: `simple` and `stitchr`.

#### `simple` mode

Simple mode outputs a single date-/version-stamped FASTA file of the requested loci/regions for a given species, as specified by the command line flag options described above.

#### `stitchr` mode

Stitchr mode downloads and formats sequence data suitable for use with [Stitchr, the tool for generation of complete coding nucleotide sequences from minimal V/J/CDR3 information](https://github.com/JamieHeather/stitchr). 

This mode downloads all loci for all regions that are available for the requested species, and saves them in a folder named after the common name of that species. That directory will contain the following files:

* `data-production-date.tsv`
    * Contains information about the IMGT and script versions used to generate this data
    * Note that as of version 0.6.1 the names of these fields has changed slightly, to allow for a more standardised metadata handling in downstream scripts
* `imgt-data.fasta`
    * Contains all of the FASTA reads that were successfully downloaded for this species
* `J-region-motifs.tsv`
    * Contains automatically inferred CDR3 junction ending motifs and residues (using the process established in [the autoDCR TCR assignation tool](https://github.com/JamieHeather/autoDCR)), for use in finding the ends of junctions in `stitchr`
* `C-region-motifs.tsv`
    * Contains automatically inferred in-frame constant region peptide sequences, for use in finding the correct frame of stitched sequences
* `TR[A/B/G/D].fasta`
    * FASTA files of the individual loci's genes

It is simply run by specifying the mode and the desired species, either using common or scientific names:

```
IMGTgeneDL -s human -n -m stitchr
IMGTgeneDL -s Homo+sapiens -m stitchr
```

The script will only output a FASTA file for a specific locus if there is sufficient information for `stitchr` to use, i.e. there must be at least one leader, variable, joining, and constant region sequence, which are not currently available for all sequences that are listed in the database. In order to be used by `stitchr` it also appends an additional field to the end of the IMGT-provided FASTA header after a '~' character, labeling it's sequence type (LEADER/VARIABLE/JOINING/CONSTANT), allowing for explicit type declaration in the face of potentially variable IMGT-provided fields.

However note that while this script will generate files if those conditions are met, these data may not be sufficient for functional `stitchr` operation without manually adding sequences (e.g. if no matching leader/variable regions are found). It's recommended that users take care when using these tools for species with relatively little banked data. 

Also note that `stitchr` mode will filter out any FASTA reads containing ambiguous or non-DNA residues, which are present for some species. All TRDV genes are also included in the `TRA.fasta` files by default, as (at least in humans) all can be found rearranged with TRAJ genes.

## Notes

### Constant region sequence downloading

The architecture of the TCR loci differs a little between the genes and across species, which the IMGT nomenclature has specific terms to cope with. However the URL based searching (or FASTA filtration in `stitchr` mode) in this script requires provision exact **exon names** for each species. However this presents an issue for constant regions, as these vart.

By default, `IMGTgeneDL` assumes the following exon configurations for each locus:

* **TRAC**: EX1+EX2+EX3+EX4UTR
* **TRBC**: EX1+EX2+EX3+EX4
* **TRGC**: EX1+EX2+EX3
* **TRDC**: EX1+EX2+EX3+EX4UTR

However requesting these complete exon arrangements automatically from IMGT often fails to produce a complete pre-spliced sequence, even when the correctly annotated individual exons appear to be present in the database. As such `IMGTgeneDL` instead downloads each exon individually and concatenates the sequences into complete constant regions as per the relevant exon orders. (This explains why '?' characters appear in the FASTA headers for constant region sequences, as this script doesn't re-calculate the relevant fields for assembled sequences, e.g. start/stop position in germline references.)

Alternative exon arrangements, as are common in gamma chain constant regions, can be specified in the `c-region-variant-configs.tsv` file. Note that in order to distinguish between potentially valid complete sequences from the same gene/allele combination but using a different arrangement of exons, `IMGTgeneDL` will append the non-standard characters to the end of the allele after an underscore. E.g. for humans, `TRGC2*05` will specify the default exon arrangement 'EX1+EX2+EX3', while `TRGC2*05_TR` will specify the sequence for 'EX1+EX2**T**+EX2**R**+EX2+EX3'. It is recommended that users intending to use `stitchr` for a locus in a species with such variant configurations use the `-p / preferred_allele_path` option in `stitchr` to specify precisely which exon configuration they need for that constant region, as neither script is capable of determining the appropriate one to use automatically.

Users should also take caution as sometimes IMGT gene annotations can err, resulting in sequences that will not combine together downstream (e.g. in `stitchr`). The most common error I've noticed is that sometimes the 3' terminal residue of the J gene which gets spliced onto the first codon of the constant region can appear in the constant region sequence inappropriately, leading to a potential frameshift. It is recommended that users inspect the `C-region-motifs.tsv` file after running this script, and pay attention to the final column: this field specifies an amino acid motif downstream of any in-frame stop codons in a constant region. Ordinarily, the only loci which should have a sequence here are TRAC and TRDC (those with a EX4UTR exon), and thus a TRBC or TRGC gene with an entry here can be suggestive of an annotation error. 

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

#### Warnings file

`IMGTgeneDL` will produce a text file whenever it is run, either named or suffixed `IMGTgeneDLwarnings.txt`. This file contains the version information, terminal command used, and any warnings generated throughout the process. IMGT/GENE-DB release number and access date are used to prefix the data (unless running in `stitchr` mode, in which case this information is provided in the `data-production-date.tsv` file). It is recommended to check these files after running, and then keep these them with the output data to maintain data provenance and provide reproducibility information when publishing. 

##### Licensing notes

The code and documentation of this repo are provided under [an MIT license](https://raw.githubusercontent.com/JamieHeather/IMGTgeneDL/main/LICENSE). I am not affiliated with IMGT, and this tool is only shared to help the field to better make of use of IMGT data. As stated [on their website](https://www.imgt.org/about/termsofuse.php) and [in their publications](https://doi.org/10.1093/nar/gkab1136), IMGT's policy regarding sharing of their data is that:

"*... IMGT® software and data are provided to the academic users and NPO's (Not for Profit Organization(s)) under the [CC BY-NC-ND 4.0 license](https://creativecommons.org/licenses/by-nc-nd/4.0/). Any other use of IMGT® material, from the private sector, needs a financial arrangement with CNRS.*"

Please TCR responsibly.