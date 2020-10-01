# NewTargets

## Description

Bioinformatic sequence recovery for universal sequence capture bait kits can be substantially improved by appropriate tailoring of target files to the group under study. To enable the best possible locus recovery from Angiosperms353 capture data, we have developed an expanded target file (`mega353.fasta`) incorporating sequences from over 500 transcriptomes from the [1KP][9] project. To ensure computational efficiency, we provide the script `filter_mega353.py` to tailor the file to the needs of a specific dataset.  `filter_mega353.py` can be used to subsample the `mega353.fasta` file based on user-selected taxa or taxon groups, as defined by unique 1KP transcriptome codes, families, orders, species, or broader groups (e.g. basal eudicots). In addition, we  provide the script `BYO_transcriptomes.py`, which can be used to incorporate any transcriptome dataset in to any target file containing protein-coding sequences. These tailored and customised target files can be used directly in target-capture pipelines such as [HybPiper][8]. 

**Data files**
- `mega353.fasta` A target file for use with target enrichment datasets captured using the Angiosperms353 bait kit. 
- `filtering_options.csv` A comma-separated values file listing the options available for filtering the `mega353.fasta` file. This reference file can also be produced by the `filter_mega353.py` script (see below).

**Scripts**
- `filter_mega353.py` A script to filter the `mega353.fasta` target file.
- `BYO_transcriptomes.py` A script to add  sequences from any transcriptome dataset to any target file.
## Dependencies

Dependencies for `filter_mega353.py`
- Python 3.7 or higher
- [BioPython][4] 1.76 or higher

Dependencies for `BYO_transcriptomes.py`
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

**Quick usage:**
```
python filter_mega353.py [-h] [-filtered_target_file FILTERED_TARGET_FILE]
                         [-report_filename REPORT_FILENAME]
                         [-list_filtering_options]
                         mega353_file select_file
```

Please see the Wiki page for [filter_mega353][7] for further details.
***

### BYO_transcriptomes.py

**Quick usage:**
```
Code Block
```

Please see the Wiki page for [BYO_transcriptomes][6] for further details.

[1]: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate/ "Link to EXONERATE download page"
[2]: http://hmmer.org/ "Link to HMMER download page"
[3]: https://mafft.cbrc.jp/alignment/software/ "Link to MAFFT download page"
[4]: https://biopython.org/wiki/Download "Link to BioPython download page"
[5]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/Installing-dependencies/ "Link to Installing dependencies Wiki page"
[6]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/Adding-transcriptome-sequences-to-a-target-file-with-BYO_transcriptomes.py/ "Link to BYO_transcriptomes Wiki page"
[7]: https://github.com/chrisjackson-pellicle/NewTargets/wiki/Filtering-the-mega353.fasta-target-file/ "Link to filter_mega353 Wiki page"
[8]: https://github.com/mossmatters/HybPiper/ "Link to the HybPiper GitHub repository"
[9]: https://sites.google.com/a/ualberta.ca/onekp/ "Link to the 1000 Plants website"



