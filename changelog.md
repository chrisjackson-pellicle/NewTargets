31 March 2023

- Added a `--version` / `-v` parameter to `BYO_transcriptome.py`. The current version is 1.0.1.
- HMM searches are now performed using the command `nhhmer` rather than `hmmsearch`. `nhhmer` checks both strands by default; where a transcriptome hit is on the negative strand, the corresponding transcript is reverse complemented before alignment. 
- Added the output folder `00_transcriptomes_renamed`; the renamed transcriptome files are now writen to this folder, rather than being written to the input transcriptome directory.
- Added option `--symfrac 0.0` to `hmmbuild` command. 


16 May 2022

- Updated the filter_mega353.py script filter list to include all `Target_name` IDs present in the current
  `mega353.fasta` file. Replaced the `filtering_options.csv` file containing 566 `Target_name` rows with the updated 
  list containing 707 `Target_name` rows.