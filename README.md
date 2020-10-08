# NewTargets

## Description

Bioinformatic sequence recovery for universal target-capture bait kits can be [substantially improved][12] by appropriate tailoring of target files to the group under study. To enable the best possible locus recovery from [Angiosperms353][10] capture data, we have developed an expanded target file (`mega353.fasta`) incorporating sequences from over 550 transcriptomes from the [1KP][9] project. To maximise computational efficiency we provide the script `filter_mega353.py`, which can be used to subsample the `mega353.fasta` file based on user-selected taxa or taxon groups, as defined by unique 1KP transcriptome codes, species, families, orders, or broader groups (e.g. Basal Eudicots, Monocots, etc). In addition, we  provide the script `BYO_transcriptome.py`, which can be used to add sequences from any transcriptome to any protein-coding nucleotide target file. These tailored and customised target files can be used directly in target-capture pipelines such as [HybPiper][8]. 

**Data files**
- `mega353.fasta` A target file for use with target enrichment datasets captured using the Angiosperms353 bait kit. 
- `filtering_options.csv` A comma-separated values file listing the options available for filtering the `mega353.fasta` file. This reference file can also be produced by the `filter_mega353.py` script (see below).

**Scripts**
- `filter_mega353.py` A script to filter the `mega353.fasta` target file.
- `BYO_transcriptome.py` A script to add sequences from any transcriptome dataset to any target file containing protein-coding sequences.

**Manuscript** 
- https://www.biorxiv.org/content/10.1101/2020.10.04.325571v1

## Dependencies

Dependencies for `filter_mega353.py`
- Python 3.7 or higher
- [BioPython][4] 1.76 or higher
- [pandas][11] 1.0.3 or higher

Dependencies for `BYO_transcriptome.py`
- Python 3.7 or higher
- [EXONERATE][1] 2.4.0
- [HMMER][2] 3.2.1 or higher
- [MAFFT][3] 7.407 or higher
- [BioPython][4] 1.76 or higher

Please see the Wiki page [Installing dependencies][5] for further details.

## Installation

Assuming all dependencies are installed, either:

1. Download the NewTargets package directly from the repository home page and unzip it. Note that the `mega353.fasta` file in the unzipped package is also provided as a `.zip` file, and will need to be unzipped separately. 
2. Clone the repository using the command `git clone https://github.com/chrisjackson-pellicle/NewTargets.git`. Unzip the `mega353.zip` file.


## Scripts

### filter_mega353.py

*Input*:
- `mega353.fasta`. The expanded Angiosperms353 target file.
- `select_file`. A text file containing a list of IDs; sequences from these samples will be retained in the filtered output target file.   

*Output*:
- `filtered_target_file`. The filtered target file containing the default Angiosperms353 sequences and any additional sequences corresponding to IDs in the `select_file`.
- `report_file`. A report file in `.csv` format, listing samples with sequences retained in the filtered target file (excluding default Angiosperms353 samples).

**Quick usage:**
```
python filter_mega353.py [-h] [-filtered_target_file FILTERED_TARGET_FILE]
                         [-report_filename REPORT_FILENAME]
                         [-list_filtering_options]
                         mega353_file select_file
```
Example command line:

`python filter_mega353.py mega353.fasta select_asparagales.txt -filtered_target_file asparagales_targetfile.fasta -report_filename asparagales_report.csv`

To generate the `filtering_options.csv` file:

`python filter_mega353.py -list_filtering_options`

Please see the Wiki page [filter_mega353][7] for further details.
***

### BYO_transcriptome.py

*Input*:
- `target_file`. A target file containing protein-coding nucleotide sequences in fasta format.
- `transcriptomes_folder`. A directory containing one or more transcriptomes in fasta format. 

*Output*:

These are the main results and reports folders; see the Wiki page for full output details.

- `17_mega_target_file`. A folder containing the final target file `BYO_target.fasta`.
- `18_reports`. A folder containing the general report file `summary_report.csv`, and the file `report_per_gene.csv` containing a presence/absence matrix of transcriptome hits for each gene/transcriptome.

**Quick usage:**
```
python BYO_transcriptome.py [-h] [-num_hits_to_recover <integer>]
                            [-python_threads <integer>]
                            [-external_program_threads <integer>]
                            [-length_percentage <float>]
                            [-hmmsearch_evalue <number in scientific notation; default is 1e-50>]
                            [-trim_to_refs <taxon_name>] [-no_n]
                            [-skip_exonerate_frameshift_fix]
                            target_file transcriptomes_folder
```
Example command line:

`python BYO_transcriptome.py asparagales_targetfile.fasta additional_asparagales_transcriptomes_folder -python_threads 4 -external_program_threads 4`

Please see the Wiki page [BYO_transcriptome][6] for further details.

[1]: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate/ "Link to EXONERATE download page"
[2]: http://hmmer.org/ "Link to HMMER download page"
[3]: https://mafft.cbrc.jp/alignment/software/ "Link to MAFFT download page"
[4]: https://biopython.org/wiki/Download "Link to BioPython download page"
[5]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/Installing-script-dependencies "Link to Installing dependencies Wiki page"
[6]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/BYO_transcriptome.py:-adding-transcriptome-sequences-to-a-target-file "Link to BYO_transcriptome Wiki page"
[7]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/filter_mega353.py:-filtering-the-mega353.fasta-target-file "Link to filter_mega353 Wiki page"
[8]: https://github.com/mossmatters/HybPiper/ "Link to the HybPiper GitHub repository"
[9]: https://sites.google.com/a/ualberta.ca/onekp/ "Link to the 1000 Plants website"
[10]: https://dx.doi.org/10.1093%2Fsysbio%2Fsyy086 "Link to Angiosperms353 manuscript"
[11]: https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html "Link to pandas installation instructions"
[12]: https://www.biorxiv.org/content/10.1101/2020.10.04.325571v1 "Link to NewTargets manuscript"



