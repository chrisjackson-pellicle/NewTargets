# NewTargets

## Description

Here we provide...

**Data files**
- `mega353.fasta` A target file for use with target enrichment data that has been captured using the Angiosperms353 bait kit. 
- `filtering_options.csv` A `.csv` file listing the options available for filtering the `mega353.fasta` file. 

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

1. Download the NewTargets package directly from the repository home page and unzip it. Note that the `mega353.fasta` file is provided as a `.zip` file, and will need to be unzipped separately. 
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

***

Please see the Wiki page for [filter_mega353][7] for further details.

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



