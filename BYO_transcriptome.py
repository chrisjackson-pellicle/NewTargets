#!/usr/bin/env python

# Author: Chris Jackson 2020

"""
Script to add sequences to the gene references in a nucleotide bait-capture target file, using a user-provided
set of transcriptomes.

########################################################################################################################

Additional information:

NOTE:

1)  Your target genes should be grouped and differentiated by a suffix in the fasta header,
    consisting of a dash followed by an ID unique to each gene, e.g.:

    >AJFN-4471
    >Ambtr-4471
    >Ambtr-4527
    >Arath-4691
    >BHYC-4691
    etc...

2)  Your transcriptomes should be named with a unique identifier; this can optionally be followed by a dash and
    additional text,
    e.g.:

    AVJK.fa.gz
    PPPZ.fa.gz
    etc...

    or

    AVJK-SOAPdenovo-Trans-assembly.fa.gz
    PPPZ-SOAPdenovo-Trans-assembly.fa.gz
    etc...

3)  In the summary report, the sum of new sequences added might differ slightly from the difference between the
    total number of sequences in the original vs new target files, if the transcriptomes searched already have sequences
    present in the target file provided.

4)  If you are providing your own target taxon names for use in transcriptome hit trimming and frameshift correction
    (using the -trim_to_refs <taxon_name> flag), you need to provide each taxon name separately e.g.:

    -trim_to_refs species1_name -trim_to_refs species2_name

########################################################################################################################

"""

import sys
try:
    import Bio
except ImportError:
    sys.exit(f"Required Python package 'Bio' not found. Is it installed for the Python used to run this script?")
try:
    import numpy
except ImportError:
    sys.exit(f"Required Python package 'numpy' not found. Is it installed for the Python used to run this script?")
import logging
import fnmatch
import sys
import argparse
import os
import socket
from Bio.Align.Applications import MafftCommandline
import subprocess
import re
import shutil
import textwrap
import itertools
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import wait
from Bio import SeqIO, AlignIO, SearchIO
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from collections import defaultdict, OrderedDict
import glob
import gzip
import bz2
import zipfile
import numpy as np
from operator import itemgetter
from multiprocessing import Manager

########################################################################################################################
########################################################################################################################
# Get current working directory and host name

cwd = os.getcwd()
host = socket.gethostname()

# Configure logger:

# Create a custom logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Create handlers
c_handler = logging.StreamHandler(sys.stdout)
existing_log_file_numbers = [int(file.split('_')[-1]) for file in os.listdir('.') if
                             fnmatch.fnmatch(file, '*.mylog*')]
if not existing_log_file_numbers:
    new_log_number = 1
else:
    new_log_number = sorted(existing_log_file_numbers)[-1] + 1
f_handler = logging.FileHandler(f'BYO_transcriptomes.mylog_{new_log_number}', mode='w')

# Create formatters and add them to handlers
c_format = logging.Formatter('%(message)s')
f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
f_handler.setFormatter(f_format)

# Add handlers to the logger
logger.addHandler(c_handler)
logger.addHandler(f_handler)

########################################################################################################################
########################################################################################################################
# Define functions:


def createfolder(directory):
    """
    Attempts to create a directory with the name provided if it doesn't exist, and provides an error message on
    failure.
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        logger.info(f'Error: Creating directory: {directory}')


def file_exists_and_not_empty(file_name):
    """
    Checks if file exists and is not empty by confirming that its size is not 0 bytes.
    """
    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def unzip(file):
    """
    Unzips a .zip file unless unzipped file already exists.
    """
    expected_unzipped_file = re.sub('.zip', '', file)
    directory_to_unzip = os.path.dirname((file))
    if not file_exists_and_not_empty(expected_unzipped_file):
        with zipfile.ZipFile(file) as infile:
            infile.extractall(directory_to_unzip)
        os.remove(file)


def gunzip(file):
    """
    Unzips a .gz file unless unzipped file already exists.
    """
    expected_unzipped_file = re.sub('.gz', '', file)
    if not file_exists_and_not_empty(expected_unzipped_file):
        with open(expected_unzipped_file, 'w') as outfile:
            with gzip.open(file, 'rt') as infile:
                outfile.write(infile.read())
        os.remove(file)


def decompress_bz2(file):
    """
    Unzips a .bz2 file unless unzipped file already exists.
    """
    expected_unzipped_file = re.sub('.bz2', '', file)
    if not file_exists_and_not_empty(expected_unzipped_file):
        with open(expected_unzipped_file, 'wb') as outfile:
            with bz2.BZ2File(file, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
        os.remove(file)


def concatenate_small(outfile, *args):
    """
    Takes an output filename and one or more files as parameters; concatenates files. Note that this will append if
    the output file already exists.
    e.g. concatenate_small('test.fastq', 'IDX01_S1_L001_R1_001.fastq', 'IDX01_S1_L001_R2_001.fastq').
    """
    with open(outfile, 'a+') as outfile:
        for filename in args:
            with open(filename, 'r') as infile:
                outfile.write(infile.read())


def pad_seq(sequence):
    """
    Pads a sequence Seq object to a multiple of 3 with 'N'.
    """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


def done_callback(future_returned):
    """
    Callback function for ProcessPoolExecutor futures; gets called when a future is cancelled or 'done'.
    """
    if future_returned.cancelled():
        logger.info(f'{future_returned}: cancelled')
    elif future_returned.done():
        error = future_returned.exception()
        if error:
            logger.info(f'\n{future_returned}: error returned: {error}')
            return
        else:
            result = future_returned.result()
    return result


def check_dependencies(target_file, transcriptomes_folder, python_threads, external_program_threads, length_percentage,
                       hmmsearch_evalue, refs_for_manual_trimming, no_n):
    """
    Checks Python version, successful import of third party Python modules, the presence of external executables, and
    HMMER, Exonerate and MAFFT versions. Prints and logs the run parameters.
    """

    # Check python version
    if sys.version_info[0:2] < (3, 7):
        sys.exit(
            f'Must be using Python 3.7 or higher. You are using version {sys.version_info[0]}.{sys.version_info[1]}. '
            f'You could try running the script in a conda environment with version 3.7 installed...')

    # Check external executables
    executables = ['hmmbuild', 'hmmsearch', 'mafft', 'exonerate']
    logger.info('')
    for executable in executables:
        if not shutil.which(executable):
            sys.exit(f"Required executable '{executable}' not found. Is it installed and in your PATH?")
        else:
            logger.info(f"Found required executable: {executable:<15}{shutil.which(executable):<30}")
    logger.info('')

    # Check HMMer version is at least 3.2.1
    hmmmer_help = subprocess.run(['hmmsearch', '-h'], capture_output=True)
    version = re.search(r'\bHMMER [0-9][.][0-9].[0-9]\b|\bHMMER [0-9][.][0-9]\b', str(hmmmer_help)).group(0)
    if not version:
        logger.info(f"Can't determine HMMer version. Please make sure you are using at least v3.2.1!\n")
    elif [int(integer) for integer in re.split(r'[.]| ', version)[1:]] < [3, 2, 1]:
        sys.exit(f'Please use at least HMMer version 3.2.1. You are using {version}')
    else:
        logger.info(f"You're using HMMer version '{version}', continuing...")

    # Check Exonerate version is at least 2.4.0
    exonerate_help = subprocess.run(['exonerate', '-h'], capture_output=True)
    version = re.search(r'\bexonerate version [0-9][.][0-9].[0-9]\b', str(exonerate_help)).group(0)
    if not version:
        logger.info(f"Can't determine Exonerate version. Please make sure you are using at least v2.4.0!\n")
    elif [int(integer) for integer in re.split(r'[.]| ', version)[2:]] < [2, 4, 0]:
        sys.exit(f'Please use at least Exonerate version 2.4.0. You are using {version}')
    else:
        logger.info(f"You're using Exonerate version '{version}', continuing...")

    # Check mafft version is at least 7.407
    mafft_version = subprocess.run(['mafft', '--version'], capture_output=True)
    version = re.search(r'\bv[0-9][.][0-9][0-9][0-9]\b', str(mafft_version)).group(0)
    if not version:
        logger.info(f"Can't determine mafft version. Please make sure you are using at least v7.407!\n")
    elif [int(integer) for integer in re.split(r'[v.]| ', version)[1:]] < [7, 407]:
        sys.exit(f'Please use at least mafft version 7.450. You are using {version}')
    else:
        logger.info(f"You're using mafft version '{version}', continuing...\n")

    # Check python modules
    python_modules = ['Bio', 'numpy']
    for module in python_modules:
        if module in sys.modules:
            logger.info(f"Imported required python module: {module:8}")

    # Check hmmsearch eValue is in scientific notation
    if not re.search(r'\b[0-9]+e-[0-9]+\b', hmmsearch_evalue):
        sys.exit(f'The eValue for hmmsearch is not in scientific notation. Value used: {hmmsearch_evalue}')

    # Print parameters for run
    default_refs = ['Ambtr', 'Arath', 'Orysa']
    try:
        commandline_refs_for_manual_trimming = list(sorted(set(refs_for_manual_trimming)))
    except TypeError:
        commandline_refs_for_manual_trimming = None
    if commandline_refs_for_manual_trimming:
        refs_for_manual_trimming = commandline_refs_for_manual_trimming
    else:
        refs_for_manual_trimming = default_refs

    logger.info(f'\n{"*"*28} Running analysis with the following settings: {"*"*29}\n\n'
                f'{"Target file:":<50}{target_file:<50}\n'
                f'{"Transcriptomes folder:":<50}{transcriptomes_folder:<50}\n'
                f'{"Python multiprocessing pool threads:":<50}{python_threads:<50}\n'
                f'{"External program threads:":<50}{external_program_threads:<50}\n'
                f'{"Length percentage cut-off for grafting:":<50}{length_percentage:<50}\n'
                f'{"eValue for hmmsearch:":<50}{hmmsearch_evalue:<50}\n'
                f'{"References for trimming transcriptomes hits:":<50}{", ".join(refs_for_manual_trimming):<50}\n'
                f'{"Remove all n characters from transcriptome hits:":<50}{str(no_n):<50}\n')
    return


def check_files_for_processing(target_fasta_file, transcriptomes_folder, refs_for_manual_trimming):
    """
    Checks the number, type and naming conventions of files and target-file fasta headers. Unzips transcriptome files
    if necessary and creates a copy of each with transcripts renamed to 'contig_1-transcriptomeID,
    contig_2-transcriptomeID' etc. Checks that a reference sequence (from either the default list: ['Ambtr', 'Arath',
    'Orysa'] or as provided by the -trim_to_refs command line option) is available for each gene in the target file.
    Checks that each sequence in the target file can be translated without stop codons (or with only a single stop
    codon found within the last 10 amino-acids) in one of the forwards frames.
    """
    transcriptomes_folder_base = os.path.basename(transcriptomes_folder)
    target_fasta_file_base = os.path.basename(target_fasta_file)

    # Check target-file fasta header formatting:
    gene_lists = defaultdict(list)
    with open(target_fasta_file, 'r') as target_file:
        seqs = SeqIO.parse(target_file, 'fasta')
        incorrectly_formatted_fasta_headers = []
        for seq in seqs:
            if not re.match('.+-[^-]+', seq.name):
                incorrectly_formatted_fasta_headers.append(seq.name)
            gene_id = re.split('-', seq.name)[-1]
            gene_lists[gene_id].append(seq)
    if incorrectly_formatted_fasta_headers:
        sys.exit(f'The target-file provided "{target_fasta_file_base}" contains the following genes with incorrectly '
                 f'formatted fasta headers: {", ".join(incorrectly_formatted_fasta_headers)}.')

    gene_names_in_target_file = [seq.name for gene_seq_list in gene_lists.values() for seq in gene_seq_list]

    # Check that seqs in target_fasta_file can be translated in one of the forwards frames:
    seqs_with_frameshifts_dict = defaultdict(list)
    sequences = list(SeqIO.parse(target_fasta_file, 'fasta'))
    for sequence in sequences:
        gene_name = sequence.name.split('-')[-1]
        num_stop_codons = pad_seq(sequence.seq).translate().count('*')
        if num_stop_codons == 0:
            logger.debug(f'Translated sequence {sequence.name} does not contain any stop codons, proceeding...')
        elif num_stop_codons == 1 and re.search('[*]', str(pad_seq(sequence.seq).translate())[-10:]):
            logger.debug(f'Translated sequence {sequence.name} contains a single stop codon in the last 10 codons, '
                         f'proceeding...')
        elif num_stop_codons > 0:
            frames_with_stop_codons = 0
            for frame_start in [1, 2]:
                num_stop_codons = pad_seq(sequence[frame_start:].seq).translate().count('*')
                if not num_stop_codons:
                    logger.debug(f'Translated sequence {sequence.name} does not contain any stop codons when '
                                 f'translated in forwards frame {frame_start + 1}, proceeding...')
                    break
                elif num_stop_codons == 1 and re.search('[*]', str(pad_seq(
                        sequence[frame_start:].seq).translate())[-10:]):
                    logger.debug(f'Translated sequence {sequence.name} contains a single stop codon in the last 10 '
                                 f'codons when translated in frame {frame_start + 1}, proceeding')
                else:
                    frames_with_stop_codons += 1
            if frames_with_stop_codons == 2:  # CJJ i.e. both 2nd and 3rd frames have stop codons in them too.
                seqs_with_frameshifts_dict[gene_name].append(sequence.name)
    if len(seqs_with_frameshifts_dict) != 0:
        logger.info(f'The following target file sequences cannot be translated without stop codons in any forwards '
                    f'frame, please correct this:')
        for key, value in seqs_with_frameshifts_dict.items():
            logger.info(f'Gene {key}: {", ".join(value)}')
        sys.exit(f'Please check your target file and try again')

    # Check if a reference is present for each gene:
    try:
        refs_for_trimming_and_exonerate = list(sorted(set(refs_for_manual_trimming)))
    except TypeError:
        refs_for_trimming_and_exonerate = ['Ambtr', 'Arath', 'Orysa']

    re_compile_string = '|'.join(refs_for_trimming_and_exonerate)
    pattern = re.compile(re_compile_string)

    genes_without_refs = []
    for gene_name, gene_list in gene_lists.items():
        gene_names = [seq.name for seq in gene_list]
        if not re.search(pattern, ','.join(gene_names)):
            genes_without_refs.append(gene_name)
    if not genes_without_refs:
        logger.info(f'Target-file fasta headers look good, continuing...')
    else:
        sys.exit(f"No reference from list {refs_for_trimming_and_exonerate} found for gene(s): "
                 f"{', '.join(genes_without_refs)}. A reference for each gene is required for alignment trimming and "
                 f"Exonerate frameshift correction steps!")

    # Unzip transcriptome files if necessary
    for transcriptome in glob.glob(f'{transcriptomes_folder}/*'):
        transcriptome_id = os.path.basename(transcriptome)
        filename, file_extension = os.path.splitext(transcriptome)
        if file_extension == '.gz':
            logger.info(f'Unzipping transcriptome {transcriptome_id}...')
            gunzip(transcriptome)
        elif file_extension == '.bz2':
            logger.info(f'Unzipping transcriptome {transcriptome_id}...')
            decompress_bz2(transcriptome)
        elif file_extension == '.zip':
            logger.info(f'Unzipping transcriptome {transcriptome_id}...')
            unzip(transcriptome)

    # Rename transcriptome sequences
    logger.info(f'Renaming transcriptome sequences...')
    transcriptomes_to_process = []
    for transcriptome in [fn for fn in glob.glob(f'{transcriptomes_folder}/*') if not fn.endswith('renamed.fasta')]:
        transcriptome_id = os.path.basename(transcriptome)
        transcriptome_unique_code = re.split('[-.]', transcriptome_id)[0]
        filename, file_extension = os.path.splitext(transcriptome)
        transcriptome_renamed = f'{filename}_renamed.fasta'
        if file_exists_and_not_empty(transcriptome_renamed):
            os.remove(transcriptome_renamed)
        transcriptomes_to_process.append(filename)
        renamed_seqs = []
        seqs = SeqIO.parse(transcriptome, 'fasta')
        contig_num = 1
        for seq in seqs:
            seq.description = seq.name
            seq.name = f'transcript_{contig_num}-{transcriptome_unique_code}'
            seq.id = f'transcript_{contig_num}-{transcriptome_unique_code}'
            renamed_seqs.append(seq)
            contig_num += 1
        with open(transcriptome_renamed, 'w') as renamed:
            SeqIO.write(renamed_seqs, renamed, 'fasta')
    if not transcriptomes_to_process:
        sys.exit(f'The folder provided "{transcriptomes_folder_base}" does not appear to contain any files!')

    logger.info(f'\n{"*" * 50}INFO{"*" * 50}\n')
    logger.info(f'{"Number of sequences in target file:":<50}{len(gene_names_in_target_file)}')
    logger.info(f'{"Number of unique gene ids in target file:":<50}{len(gene_lists)}')
    logger.info(f'{"Number of transcriptomes to search:":<50}{len(transcriptomes_to_process)}\n')
    logger.info(f'If the numbers above do not match your expectations, please check your file and sequences names.')
    logger.info(f'\n{"*" * 104}\n')
    return


def split_targets(target_fasta_file, output_folder):
    """
    Takes the target fasta file and splits it into one fasta file per gene (grouping via the uniqueID suffix in
    the fasta header e.g.

    >AJFN-4471
    >Ambtr-4471
    ...

    """
    createfolder(output_folder)
    gene_lists = defaultdict(list)
    with open(target_fasta_file, 'r') as target_file:
        seqs = SeqIO.parse(target_file, 'fasta')
        for seq in seqs:
            gene_id = re.split('-', seq.name)[-1]
            gene_lists[gene_id].append(seq)
    for key, value in gene_lists.items():
        with open(f'{output_folder}/{key}.fasta', 'w') as outfile:
            for sequence in value:
                SeqIO.write(sequence, outfile, 'fasta')


def align_targets(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Returns filename
    of the alignment produced.
    """
    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
    except AssertionError:
        mafft_cline = (MafftCommandline(algorithm, thread=threads, input=fasta_file))
        stdout, stderr = mafft_cline()
        with open(expected_alignment_file, 'w') as alignment_file:
            alignment_file.write(stdout)
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished generating alignment for file {fasta_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return os.path.basename(expected_alignment_file)


def align_targets_multiprocessing(target_gene_folder, alignments_output_folder, algorithm='linsi', pool_threads=1,
                                  mafft_threads=2):
    """
    Generate alignments via function <align_targets> using multiprocessing.
    """
    createfolder(alignments_output_folder)
    logger.info('===> Generating target gene alignments to construct HMM profiles...')
    target_genes = [file for file in sorted(glob.glob(f'{target_gene_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(align_targets, fasta_file, algorithm, alignments_output_folder, counter, lock,
                                      num_files_to_process=len(target_genes), threads=mafft_threads)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files.\n')
    if len(target_genes) != len(alignment_list):
        sys.exit(f'Only {len(alignment_list)} alignments were generated from {len(target_genes)} fasta files, check '
                  f'for errors!')


def create_hmm_profile(alignment, output_folder, counter, lock, num_files_to_process, hmmbuild_threads=2):
    """
    Create HMM profile from alignment using HMMer. Returns filename of the HMM profile produced.
    """
    createfolder(output_folder)
    alignment_file_basename = os.path.basename(alignment)
    expected_hmm_file = f'{output_folder}/{re.sub(".aln.fasta", ".hmm", alignment_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_hmm_file)
        logger.debug(f'HMM profile exists for {alignment_file_basename}, skipping...')
    except AssertionError:
        subprocess.run(['hmmbuild', '-o', '/dev/null', '--cpu', str(hmmbuild_threads), '--dna', expected_hmm_file,
                        alignment], check=True)
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished generating HMM profile for alignment {alignment_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return alignment_file_basename


def create_hmm_profile_multiprocessing(alignments_folder, hmm_output_folder, pool_threads=1, hmmbuild_threads=2):
    """
    Generate HMM profiles via function <create_hmm_profile> using multiprocessing.
    """
    createfolder(hmm_output_folder)
    logger.info('===> Generating hmm profiles from gene alignments...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(create_hmm_profile, alignment, hmm_output_folder, counter, lock,
                                      num_files_to_process=len(alignments), hmmbuild_threads=hmmbuild_threads)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    hmm_list = [hmm for hmm in glob.glob(f'{hmm_output_folder}/*.hmm') if file_exists_and_not_empty(hmm)]
    logger.info(f'\n{len(hmm_list)} HMM profiles generated from {len(future_results)} alignment files.\n')
    if len(alignments) != len(hmm_list):
        sys.exit(f'Only {len(hmm_list)} hmm profiles were generated from {len(alignments)} alignments, check '
                 f'for errors!')


def hmm_vs_transcriptome(hmm_profile, transcriptomes_folder, hits_output_folder, hmm_logs_output_folder,
                         num_hits_to_recover, counter, lock, num_files_to_process, hmmsearch_threads=1,
                         hmmsearch_evalue='1e-50'):
    """
    Performs HMM searches against transcriptome fasta files using the provided HMM profile, extracts hits from
    transcriptome and writes them to a fasta file.
    """
    hmm_profile_basename = os.path.basename(hmm_profile)
    gene_name = hmm_profile_basename.split('.')[0]
    logger.debug(f'SEARCH: searching with profile {hmm_profile_basename}')
    createfolder(f'{hits_output_folder}/{gene_name}')

    # If a HMM hit file doesn't already exist, search each transcriptome with the HMM profile
    for transcriptome in sorted(glob.glob(f'{transcriptomes_folder}/*renamed.fasta'),
                                key=lambda a: re.split('/', a)[-1]):
        transcriptome_id = os.path.basename(transcriptome)
        fasta_hit_file = f'{hits_output_folder}/{gene_name}/{transcriptome_id}_{gene_name}_HMM_hits.fasta'
        seqs_to_recover = []
        if file_exists_and_not_empty(fasta_hit_file):
            logger.debug(f'A fasta file of HMM hits for gene {gene_name} vs transcriptome {transcriptome_id} already '
                         f'exists, skipping...')
            continue
        logger.debug(f'SEARCH: searching transcriptome {transcriptome_id} with profile {hmm_profile_basename}')
        log_file = f'{hmm_logs_output_folder}/{hmm_profile_basename}_{transcriptome_id}.log'
        human_readable_log_file = f'{hmm_logs_output_folder}/' \
                                  f'{hmm_profile_basename}_{transcriptome_id}.human_readable.log'
        subprocess.run(['hmmsearch', '--tblout', log_file, '-E', hmmsearch_evalue, '-o', human_readable_log_file,
                        '--cpu', str(hmmsearch_threads), hmm_profile, transcriptome], check=True)

        # Read in log file and recover hit name(s)
        with open(log_file, 'r') as hmm_results:
            results = hmm_results.readlines()
            hits = [line.split()[0] for line in results if not line.startswith('#')]
        if hits:
            logger.debug(f'*** Transcriptome {transcriptome_id} contains {len(hits)} hits for '
                         f'{hmm_profile_basename} ***')
            for hit in hits[0:num_hits_to_recover]:
                seqs_to_recover.append(hit)
        else:
            logger.debug(f'No hits for {gene_name} from transcriptome {transcriptome_id}...')
            pass

        # Parse transcriptome fasta file and recover hit sequences
        with open(transcriptome) as transcriptome_fasta:
            seq_records_to_write = []
            seqs = SeqIO.parse(transcriptome_fasta, 'fasta')
            for seq in seqs:
                if seq.name in seqs_to_recover:
                    seq_records_to_write.append(seq)

        # Write the hits sequences to a fasta file
        with open(fasta_hit_file, 'w') as fasta_output:
            for seq_record in seq_records_to_write:
                SeqIO.write(seq_record, fasta_output, 'fasta')
    with lock:
        counter.value += 1
        print(f'\rFinished searching transcriptomes with {hmm_profile_basename}, '
              f'{counter.value}/{num_files_to_process}', end='')
    return hmm_profile_basename


def hmm_vs_transcriptome_multiprocessing(hmmprofile_folder, transcriptomes_folder, hits_output_folder,
                                         hmm_logs_output_folder, num_hits_to_recover, pool_threads=1,
                                         hmmsearch_threads=1, hmmsearch_evalue='1e-50'):
    """
    Run HMM searches of transcriptomes via function <hmm_vs_transcriptome> using multiprocessing.
    """
    createfolder(hits_output_folder)
    createfolder(hmm_logs_output_folder)
    logger.info('===> Searching transcriptomes with HMM profiles...')

    hmm_profiles = [hmm for hmm in sorted(glob.glob(f'{hmmprofile_folder}/*.hmm'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(hmm_vs_transcriptome, hmm_profile, transcriptomes_folder, hits_output_folder,
                                      hmm_logs_output_folder, num_hits_to_recover, counter, lock,
                                      num_files_to_process=len(hmm_profiles), hmmsearch_threads=hmmsearch_threads,
                                      hmmsearch_evalue=hmmsearch_evalue)
                          for hmm_profile in hmm_profiles]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")


def align_extractions_single_reference(single_gene_alignment, single_gene_alignment_object, concatenated_hits,
                                       mafft_threads, single_gene_alignment_with_hits_name):
    """
    Processes alignments of transcriptome hits to a target file reference in the case where only a single reference is
    present in the target file.
    """
    for single_sequence in single_gene_alignment_object:
        single_sequence.id = f'_seed_{single_sequence.name}'
        single_sequence.description = ''
        single_sequence.name = ''
        single_ref_sequence = single_sequence

    if file_exists_and_not_empty(concatenated_hits):
        single_gene_alignment = re.sub('.fasta', '.single_seed.fasta', single_gene_alignment)
        seqs_to_align = [single_ref_sequence]
        for hit_seq in SeqIO.parse(concatenated_hits, 'fasta'):
            seqs_to_align.append(hit_seq)
        with open(single_gene_alignment, 'w') as seed_single_file:
            SeqIO.write(seqs_to_align, seed_single_file, 'fasta')
        mafft_cline = (MafftCommandline(localpair=True, thread=mafft_threads, input=single_gene_alignment,
                                        lop=-5.00, lep=-0.5, lexp=-0.1))
        stdout, stderr = mafft_cline()
        with open(single_gene_alignment_with_hits_name, 'w') as alignment_file:
            alignment_file.write(stdout)
        os.remove(single_gene_alignment)

    elif not file_exists_and_not_empty(concatenated_hits):  # i.e. if there are no transcriptome hits.
        with open(single_gene_alignment_with_hits_name, 'w') as alignment_file:
            SeqIO.write(single_ref_sequence, alignment_file, 'fasta')


def strip_n(concatenated_hits):
    """
    Strips the character N from a fasta file of nucleotide sequences.
    """
    concatenated_seqs = SeqIO.parse(concatenated_hits, 'fasta')
    stripped_seqs_to_write = []
    for seq in concatenated_seqs:
        seq.seq = seq.seq.ungap(gap='N')
        stripped_seqs_to_write.append(seq)
    with open(concatenated_hits, 'w') as n_stripped_concatenated_hits:
        SeqIO.write(stripped_seqs_to_write, n_stripped_concatenated_hits, 'fasta')


def align_extractions(single_gene_alignment, output_folder, hit_folder, concatenated_folder, seqs_with_ns_folder,
                      counter, lock, num_files_to_process, mafft_threads=2, no_n=False):
    """
    Takes a single gene alignment (from folder_02) from the target fasta file and aligns the hits extracted from
    the transcriptomes (from folder_04) using the mafft 'seed' option, which prefixes sequence fasta headers in the
    original alignment with the string '_seed_'. Note that this occurs even if there are no sequences in the
    'concatenated hits' file. In cases where the target file only contains a single sequence for a given gene, the
    string '_seed_' is manually added as a prefix in the fasta header, and a standard alignment is performed. In the
    latter case, if there are no sequences in the 'concatenated hits' file, the prefixed target sequence is written to
    to a new file by itself.
    """
    single_gene_alignment_object = AlignIO.read(single_gene_alignment, 'fasta')
    seqs = [seq for seq in single_gene_alignment_object]
    if len(seqs) == 1:
        single_reference = True
    else:
        single_reference = False

    single_gene_alignment_name = os.path.basename(single_gene_alignment)
    single_gene_alignment_with_hits_name = f'{output_folder}/' \
                                           f'{re.sub(".aln.fasta", ".aln_added.fasta", single_gene_alignment_name)}'
    gene_id = (single_gene_alignment_name.replace('.aln.fasta', ''))
    seqs_with_ns_folder_gene_folder = f'{seqs_with_ns_folder}/{gene_id}'

    # Create a dictionary of transcriptome hit sequences containing Ns, write sequences to file:
    concatenated_hits = f'{concatenated_folder}/{gene_id}.concat.fasta'  # Write new concatenated hits file
    if file_exists_and_not_empty(concatenated_hits):
        os.remove(concatenated_hits)
    for fasta_file in glob.glob(f'{hit_folder}/{gene_id}/*'):
        concatenate_small(concatenated_hits, fasta_file)
    seqs_with_n_dict = defaultdict(list)  # Create a dictionary of seqs that contain Ns
    for seq in SeqIO.parse(concatenated_hits, 'fasta'):
        n_count = seq.seq.count('N')
        if n_count:
            seqs_with_n_dict[gene_id].append(seq)
            createfolder(seqs_with_ns_folder_gene_folder)
            with open(f'{seqs_with_ns_folder_gene_folder}/{gene_id}_seqs_with_ns.fasta', 'w') as seqs_ns:
                SeqIO.write(seqs_with_n_dict[gene_id], seqs_ns, 'fasta')

    try:
        assert file_exists_and_not_empty(single_gene_alignment_with_hits_name)
        logger.debug(f' Alignment exists for {single_gene_alignment_name}, skipping...')
    except AssertionError:
        if no_n and file_exists_and_not_empty(concatenated_hits):  # If requested, strip Ns from transcriptome hits
            strip_n(concatenated_hits)
        if single_reference:  # i.e. if there's only a single target sequence for this gene in the target file.
            align_extractions_single_reference(single_gene_alignment, single_gene_alignment_object, concatenated_hits,
                                               mafft_threads, single_gene_alignment_with_hits_name)
        elif not single_reference:
            mafft_cline = (MafftCommandline(localpair=True, thread=mafft_threads, input=concatenated_hits, lop=-5.00,
                                            lep=-0.5, lexp=-0.1, seed=single_gene_alignment))
            stdout, stderr = mafft_cline()
            with open(single_gene_alignment_with_hits_name, 'w') as alignment_file:
                alignment_file.write(stdout)
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished aligning transcriptome hits for {single_gene_alignment_name}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        if seqs_with_n_dict:
            return single_gene_alignment_name, seqs_with_n_dict
        else:
            return single_gene_alignment_name


def align_extractions_multiprocessing(alignments_folder, output_folder, hit_folder, seqs_with_ns_folder,
                                      concatenated_folder, mafft_threads=2, pool_threads=1, no_n=False):
    """
    Generate alignments via function <align_extractions> using multiprocessing.
    """
    createfolder(output_folder)
    createfolder(concatenated_folder)
    logger.info('\n\n===> Adding transcriptome hits to gene alignments...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*aln.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(align_extractions, alignment, output_folder, hit_folder, concatenated_folder,
                                      seqs_with_ns_folder, counter, lock, num_files_to_process=len(alignments),
                                      mafft_threads=mafft_threads, no_n=no_n) for
                          alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

        seqs_with_ns = []
        for future in future_results:
            try:
                single_gene_alignment_name, seqs_with_n_dict = future.result()
                seqs_with_ns.append(seqs_with_n_dict)
            except ValueError:
                single_gene_alignment_name = future.result()

        alignment_list = [alignments for alignment in glob.glob(f'{output_folder}/*.aln_added.fasta') if
                          file_exists_and_not_empty(alignment)]
        logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} alignment files.\n')
        if len(alignments) != len(alignment_list):
            sys.exit(f'Only {len(alignment_list)} alignments were generated from {len(alignments)} genes, check '
                     f'for errors!')
        shutil.rmtree(concatenated_folder)
        return seqs_with_ns


def trim_alignments_manually(gene_alignment, output_folder, refs_for_trimmed):
    """
    Takes a fasta alignment file as input, trims it to the longest of sequences from default (_seed_Arath, _seed_Ambtr,
    _seed_Orysa), or a user provided list of taxon names.
    """
    try:
        commandline_refs_for_trimming = list(sorted(set(refs_for_trimmed)))
    except TypeError:
        commandline_refs_for_trimming = None
    if commandline_refs_for_trimming:
        sorted_refs_for_trimming = [f'_seed_{name}' for name in commandline_refs_for_trimming]
        logger.debug(f'Non-default references used for manual trimming...')
        re_compile_string = '|'.join(sorted_refs_for_trimming)
    else:
        re_compile_string = '_seed_Ambtr|_seed_Arath|_seed_Orysa'
        logger.debug('Default references used for manual trimming')

    alignment_name = os.path.basename(gene_alignment)
    trimmed_alignment = re.sub('.aln_added.*', '.aln.trimmed.fasta', alignment_name)
    expected_output_file = f'{output_folder}/{trimmed_alignment}'

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Trimmed alignment already exists for {alignment_name}, skipping....')
    except AssertionError:
        untrimmed_alignment = AlignIO.read(gene_alignment, "fasta")
        pattern = re.compile(re_compile_string)
        seed_names_and_lengths = [(seq.name, len(seq.seq.ungap(gap='-'))) for seq in untrimmed_alignment if
                                  re.search(pattern, seq.name)]
        longest_ref_name = max(seed_names_and_lengths, key=itemgetter(1))[0]
        reference_index_count = 0
        trimmed_ref_found = False
        for sequence in untrimmed_alignment:
            if sequence.name == longest_ref_name:
                trimmed_ref_found = True
                break
            else:
                reference_index_count += 1  # get index for reference seq to trim to
        if trimmed_ref_found:
            logger.debug(f'Found reference {longest_ref_name} at index {reference_index_count}!')
        untrimmed_array = np.array([sequence.seq for sequence in untrimmed_alignment])
        five_prime_slice = 0
        for position in untrimmed_array.T:
            if "-" in position[reference_index_count]:
                five_prime_slice += 1
            else:
                break
        three_prime_slice = 0
        for position in np.flipud(untrimmed_array.T):  # Transpose the numpy array
            if "-" in position[reference_index_count]:
                three_prime_slice += 1
            else:
                break
        corrected_three_prime_slice = untrimmed_alignment.get_alignment_length() - three_prime_slice
        sliced_alignment = untrimmed_alignment[:, five_prime_slice:corrected_three_prime_slice]
        with open(f'{output_folder}/{trimmed_alignment}', 'w') as outfile:
            AlignIO.write(sliced_alignment, outfile, 'fasta')
    finally:
        return expected_output_file


class ExtendDistanceCalculator(DistanceCalculator):
    """
    Subclass DistanceCalculator, adding a [0,1] based scoring matrix for _pairwise, allowing identity comparisons
    using self.skip_letters.
    """
    cjj = [[1], [0, 1], [0, 0, 1], [0, 0, 0, 1]]
    blastn = [[5], [-4, 5], [-4, -4, 5], [-4, -4, -4, 5]]
    trans = [[6], [-5, 6], [-5, -1, 6], [-1, -5, -5, 6]]
    dna_matrices = {"blastn": blastn, "trans": trans, "cjj": cjj}
    dna_models = list(dna_matrices.keys())
    dna_alphabet = ["a", "t", "c", "g"]


def get_graft_alignment(trimmed_dm, trimmed_alignment, sequence_to_graft):
    """
    For a given sequence, identifies the closest identity sequence in a given alignment. Returns an alignment
    of both sequences, and the name of the sequence used from grafting.
    """
    dm = trimmed_dm
    names = dm.names
    distance_values = dm[sequence_to_graft.name]
    sorted_distance_values = sorted(distance_values, key=float)
    for distance in sorted_distance_values:
        index = dm[sequence_to_graft.name].index(distance)
        name = names[index]
        if re.search('_seed', name):
            for sequence in trimmed_alignment:
                if sequence.name == name:
                    seq_to_graft_with = sequence
            break
    logger.debug(f'Sequence with highest identity to {sequence_to_graft.name} is {seq_to_graft_with.name}, grafting...')
    graft_alignment = MultipleSeqAlignment([sequence_to_graft, seq_to_graft_with])
    return graft_alignment, seq_to_graft_with.name


def write_fasta_and_mafft_align(original_alignment, trimmed_alignment, mafft_threads):
    """Realign a given alignment using mafft"""
    seqs = []
    with open(f'{trimmed_alignment}_tmp.fasta', 'w') as tmp:
        for seq in AlignIO.read(trimmed_alignment, 'fasta'):
            seq.seq = re.sub('-', '', str(seq.seq))
            seq = SeqRecord(Seq(seq.seq), id=seq.id, name=seq.name, description=seq.description)
            seqs.append(seq)
        SeqIO.write(seqs, tmp, 'fasta')
    mafft_cline = (MafftCommandline(localpair=True, lep=-1.0, thread=mafft_threads,
                                    input=f'{trimmed_alignment}_tmp.fasta'))
    stdout, stderr = mafft_cline()
    os.remove(f'{trimmed_alignment}_tmp.fasta')
    with open(original_alignment, 'w') as alignment_file:
        alignment_file.write(stdout)


def new_seqs_longer_than_seeds(alignment_object):
    """
    Returns true if the longest of the newly added sequences (i.e. those that do not have the prefix _seed_) is more
    than 1.1x the length of the longest seed sequence.
    """
    longest_seed = max(len(seq.seq.ungap(gap='-')) for seq in alignment_object if re.search('_seed', seq.name))
    try:
        longest_new_seq = max(len(seq.seq.ungap(gap='-')) for seq in alignment_object if not re.search(
            '_seed', seq.name))
        if longest_new_seq > (longest_seed * 1.1):
            return True
        else:
            return False
    except ValueError:
        return False


def trim_and_discard_or_graft(alignment, trimmed_alignment_folder, alignments_for_grafting_folder, new_sequence_folder,
                              length_percentage, counter, lock, num_files_to_process, refs_for_trimmed,
                              discard_short=False, mafft_threads=1):
    """
    Trims alignment to the longest of a given set of original target taxa (seed_Arath, seed_Ambtr or seed_
    Orysa from the 353Angiosperm target file by default, or user provided).

    Checks whether newly added sequences are less than <user_provided>% length (default 0.85) of the longest
    untrimmed original target sequence. If so, the new sequence is either discarded or grafted with the closest id
    sequence, calculated via a distance matrix. NOTE: this step is necessary for the Exonerate step of HybPiper,
    where SPAdes contigs are compared to the sequence with the highest BWA alignment score, and we don't want truncated
    sequences to be chosen.

    [Uses the functions: trim_alignments_manually, get_graft_alignment, new_seqs_longer_than_seeds,
    write_fasta_and_mafft_align]
    """
    alignment_name = os.path.basename(alignment)
    logger.debug(f'Processing alignment {alignment_name}')
    gene_name = alignment_name.split('.')[0]
    expected_output_file = f'{new_sequence_folder}/{gene_name}.fasta'
    warning = False
    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'A fasta file of new sequences to add already exists for {alignment_name}, skipping....')
        trimmed_alignment = AlignIO.read(f'{trimmed_alignment_folder}/{gene_name}.aln.trimmed.fasta', "fasta")
        if new_seqs_longer_than_seeds(trimmed_alignment):
            warning = f'***WARNING*** A newly added sequence is longer than it should be for gene {gene_name}!'
    except AssertionError:
        trimmed = trim_alignments_manually(alignment, trimmed_alignment_folder, refs_for_trimmed)

        # Read the trimmed alignment
        trimmed_alignment = AlignIO.read(trimmed, "fasta")

        # Check if any transcriptome hit sequence is more than 10% longer than the longest seed sequence and, if it
        # is, realign and re-trim
        if new_seqs_longer_than_seeds(trimmed_alignment):
            write_fasta_and_mafft_align(alignment, trimmed, mafft_threads)
            os.remove(trimmed)
            trimmed = trim_alignments_manually(alignment, trimmed_alignment_folder, refs_for_trimmed)
            trimmed_alignment = AlignIO.read(trimmed, "fasta")
            if new_seqs_longer_than_seeds(trimmed_alignment):
                warning = f'***WARNING*** A newly added sequence is longer than it should be for gene {gene_name}!'

        # Get the no-gap length of the longest of the seed sequences for this gene
        alignment_read = AlignIO.read(alignment, "fasta")
        longest_seed_seq = max(len(seq.seq.ungap(gap='-')) for seq in alignment_read if re.search('_seed', seq.name))

        # calculate a distance matrix for selecting _seed_ to graft with if necessary
        skip_letters = set(letter for sequence in trimmed_alignment for letter in sequence.seq if
                           letter not in ['a', 't', 'c', 'g'])
        my_calculator = ExtendDistanceCalculator('cjj', skip_letters=''.join(skip_letters))

        if not discard_short:
            logger.debug(f'Calculating distance matrix for trimmed alignment {alignment_name}')
            trimmed_dm = my_calculator.get_distance(trimmed_alignment)

        # Discard or graft newly added sequences if they fall beneath a percentage length cut-off
        trimmed_seqs_to_write = []
        for trimmed_seq in trimmed_alignment:
            if not re.search('_seed', trimmed_seq.name):
                trimmed_length = len(trimmed_seq.seq.ungap(gap='-'))
                length_fraction = trimmed_length / longest_seed_seq
                if discard_short:
                    if length_fraction < length_percentage:
                        logger.debug(f'Sequence {trimmed_seq.name} is less than {length_percentage} of the longest '
                                     f'target sequence for this gene, discarding...')
                    else:
                        trimmed_seqs_to_write.append(trimmed_seq)
                else:
                    if length_fraction < length_percentage:
                        new_sequence_to_graft = trimmed_seq
                        logger.debug(f'Sequence {trimmed_seq.name} is {trimmed_length} bases long, less than '
                                     f'{length_percentage} of the longest target sequence ({longest_seed_seq}) '
                                     f'for this gene. Grafting sequence with closest identity...')
                        grafted_gene_directory = f'{alignments_for_grafting_folder}/{gene_name}'
                        createfolder(grafted_gene_directory)
                        graft_alignment, seq_to_graft_name = get_graft_alignment(trimmed_dm, trimmed_alignment,
                                                                                 new_sequence_to_graft)

                        with open(f'{grafted_gene_directory}/{trimmed_seq.name}.aln.fasta', 'w') \
                                as outfile:
                            AlignIO.write(graft_alignment, outfile, 'fasta')
                        align_array = np.array([rec.seq for rec in graft_alignment])
                        grafted_seq_5prime = []
                        grafted_seq_3prime = []
                        grafted_name = f'{trimmed_seq.name}{seq_to_graft_name}'

                        for position in align_array.T:  # transpose
                            if '-' in position[0] and '-' in position[1]:  # skip over alignment if neither sequence yet
                                continue
                            elif '-' in position[0]:  # i.e. the alignment only has sequence from seq_to_graft
                                grafted_seq_5prime.append(position[1])
                            else:
                                break
                        for position in np.flipud(align_array.T):
                            if '-' in position[0] and '-' in position[1]:  # skip over alignment if neither sequence yet
                                continue
                            elif '-' in position[0]:  # i.e. the alignment only has sequence from seq_to_graft
                                grafted_seq_3prime.insert(0, position[1])
                            else:
                                break
                        grafted_seq_5prime_concat = ''.join(letter for letter in grafted_seq_5prime)
                        grafted_seq_3prime_concat = ''.join(letter for letter in grafted_seq_3prime)
                        final_grafted_seq = ''.join([grafted_seq_5prime_concat, str(trimmed_seq.seq.ungap(gap='-')),
                                                     grafted_seq_3prime_concat])
                        final_grafted_seq_obj = SeqRecord(Seq(final_grafted_seq), id=gene_name, name=grafted_name,
                                                          description="grafted sequence")
                        trimmed_seqs_to_write.append(final_grafted_seq_obj)
                    else:
                        trimmed_seqs_to_write.append(trimmed_seq)

        # If there are new sequences to add to the target file, write them to a fasta file
        if trimmed_seqs_to_write:
            seqfile_to_write = f'{new_sequence_folder}/{gene_name}.fasta'
            with open(seqfile_to_write, 'w') as seqfile:
                for sequence in trimmed_seqs_to_write:
                    sequence.seq = sequence.seq.ungap(gap='-')
                    transcriptome_id = re.split('-|_', sequence.name)[2]
                    grafted_gene_id = ''.join(re.split('_', sequence.name)[-1])
                    if sequence.description == "grafted sequence":
                        sequence.name = f'{transcriptome_id}_grafted_with_{grafted_gene_id}'
                    else:
                        sequence.name = f'{transcriptome_id}-{gene_name}'
                    seqfile.write(f'>{sequence.name}\n{sequence.seq}\n')
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished trimmming and grafting/discarding transcriptome hits for {alignment_name}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return alignment_name, warning

def trim_and_discard_or_graft_multiprocessing(alignments_folder, trimmed_alignment_folder,
                                              alignments_for_grafting_folder, new_sequence_folder, pool_threads,
                                              mafft_threads, length_percentage, refs_for_trimmed, discard_short=False):
    """
    Run trim_and_discard_or_graft for each alignment via function <trim_and_discard_or_graft> via multiprocessing.
    """
    createfolder(trimmed_alignment_folder)
    createfolder(new_sequence_folder)
    logger.info('===> Trimming alignment and grafting or discarding new sequences...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(trim_and_discard_or_graft, alignment, trimmed_alignment_folder,
                                      alignments_for_grafting_folder, new_sequence_folder, length_percentage, counter,
                                      lock, num_files_to_process=len(alignments), refs_for_trimmed=refs_for_trimmed,
                                      discard_short=discard_short, mafft_threads=mafft_threads)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

        long_seqs_warnings = []
        for future in future_results:
            try:
                alignment, warning = future.result()
                if warning:
                    long_seqs_warnings.append(warning)
            except ValueError:
                alignment = future.result()

    processed_alignments = [alignment for alignment in glob.glob(f'{trimmed_alignment_folder}/*.trimmed.fasta') if
                            file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(processed_alignments)} alignments processed from {len(future_results)} alignment files.\n')
    if len(alignments) != len(processed_alignments):
        sys.exit(f'Only {len(processed_alignments)} alignments were processed from {len(alignments)} genes, check '
                 f'for errors!')
    return long_seqs_warnings


def create_new_targets_file(gene_fasta_file, seqs_to_add_folder, seqs_to_write_folder):
    """
    Writes a fasta file (not an alignment) of the original target sequences and any sequences to be added.
    """
    gene_fasta_file_base = os.path.basename(gene_fasta_file)
    createfolder(seqs_to_write_folder)
    expected_output_file = f'{seqs_to_write_folder}/{gene_fasta_file_base}'

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'New target file already exists for {gene_fasta_file_base}, skipping....')
        return expected_output_file
    except AssertionError:
        seqs_to_write = []
        with open(gene_fasta_file, 'r') as single_gene:
            seqs = SeqIO.parse(single_gene, 'fasta')
            for seq in seqs:
                seqs_to_write.append(seq)
        if os.path.exists(f'{seqs_to_add_folder}/{gene_fasta_file_base}'):
            with open(f'{seqs_to_add_folder}/{gene_fasta_file_base}') as seqs_to_add:
                seqs = SeqIO.parse(seqs_to_add, 'fasta')
                for seq in seqs:
                    if seq.id not in [seq.id for seq in seqs_to_write]:
                        seqs_to_write.append(seq)
        with open(f'{seqs_to_write_folder}/{gene_fasta_file_base}', 'w') as to_write:
            for seq in seqs_to_write:
                SeqIO.write(seq, to_write, 'fasta')


def check_and_correct_reading_frames(single_gene_new_target_file, frameshifts_folder, uncorrected_frameshifts_folder,
                                     exonerate_results_folder, refs_for_exonerate, counter, lock, num_files_to_process,
                                     skip_fix_frameshifts_with_exonerate=False):
    """
    For a new gene target file, check that each sequence can be translated without stop codons from the first
    nucleotide. If not, correct the starting nucleotide. If stop codons are present regardless of whether the sequence
    is translated at position 1, 2 or 3, fix frameshifts using Exonerate. If parameter
    skip_fix_frameshifts_with_exonerate=True, simply remove the sequence and provide a warning.

    Uses the function <run_exonerate>.
    """

    target_file_basename = os.path.basename(single_gene_new_target_file)
    gene_name = f"{target_file_basename.split('.')[0]}"
    gene_with_frameshifts_folder = f'{frameshifts_folder}/{gene_name}'
    gene_with_uncorrected_frameshifts_folder = f'{uncorrected_frameshifts_folder}/{gene_name}'
    gene_exonerate_folder = f'{exonerate_results_folder}/{gene_name}'
    gene_with_frameshift_file = f'{gene_with_frameshifts_folder}/{gene_name}_with_frameshifts.fasta'
    gene_with_uncorrected_frameshifts_file = f'{gene_with_uncorrected_frameshifts_folder}/' \
                                             f'{gene_name}_with_uncorrected_frameshifts.fasta'
    seqs_with_frameshifts_dict = defaultdict(list)
    if file_exists_and_not_empty(gene_with_frameshift_file):
        os.remove(gene_with_frameshift_file)
    if file_exists_and_not_empty(gene_with_uncorrected_frameshifts_file):
        os.remove(gene_with_uncorrected_frameshifts_file)

    try:
        seqs_to_retain = []
        sequences = list(SeqIO.parse(single_gene_new_target_file, 'fasta'))
        open_frame_found = True
        for sequence in sequences:
            num_stop_codons = pad_seq(sequence.seq).translate().count('*')
            if num_stop_codons == 0:
                seqs_to_retain.append(sequence)
            elif num_stop_codons == 1 and re.search('[*]', str(pad_seq(sequence.seq).translate())[-10:]):
                logger.debug(f'Only one stop codon within the last 10 codons for seq for {sequence.name},\n '
                             f'{pad_seq(sequence.seq).translate()}')
                seqs_to_retain.append(sequence)
            elif num_stop_codons > 0:
                frames_with_stop_codons = 0
                for frame_start in [1, 2]:
                    num_stop_codons = pad_seq(sequence[frame_start:].seq).translate().count('*')
                    if not num_stop_codons:
                        sequence = sequence[frame_start:]
                        seqs_to_retain.append(sequence)
                        break
                    elif num_stop_codons == 1 and re.search('[*]', str(pad_seq(
                            sequence[frame_start:].seq).translate())[-10:]):
                        logger.debug(f'Only one stop codon within the last 10 codons for seq for {sequence.name},\n '
                                     f'{pad_seq(sequence.seq).translate()}')
                        sequence = sequence[frame_start:]
                        seqs_to_retain.append(sequence)
                    else:
                        frames_with_stop_codons += 1
                if frames_with_stop_codons == 2:  # CJJ i.e. both 2nd and 3rd frames have stop codons in them too.
                    open_frame_found = False
                    seqs_with_frameshifts_dict[gene_name].append(sequence)
        with open(single_gene_new_target_file, 'w') as checked_target_file:
            SeqIO.write(seqs_to_retain, checked_target_file, 'fasta')
        if not open_frame_found:
            warning = f'***WARNING*** Target file {target_file_basename} contains at least one sequence with no ' \
                      f'full-length ORF in any reading frame; such sequences have been removed!'
            createfolder(gene_with_frameshifts_folder)
            with open(gene_with_frameshift_file, 'w') as frameshifts:
                SeqIO.write(seqs_with_frameshifts_dict[gene_name], frameshifts, 'fasta')

            # Fix frameshifts remaining using Exonerate:
            if not skip_fix_frameshifts_with_exonerate:
                seqs_with_frameshifts_dict = run_exonerate(seqs_with_frameshifts_dict, refs_for_exonerate,
                                                           single_gene_new_target_file, gene_name,
                                                           gene_exonerate_folder)
            if seqs_with_frameshifts_dict[gene_name]:  # i.e. some seqs couldn't be fixed...
                createfolder(gene_with_uncorrected_frameshifts_folder)
                with open(gene_with_uncorrected_frameshifts_file, 'w') as uncorrected_frameshifts:
                    SeqIO.write(seqs_with_frameshifts_dict[gene_name], uncorrected_frameshifts, 'fasta')
    except:
        logger.debug(f'Checking and correcting read frames failed for gene {gene_name}')
        raise
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished checking and correcting reading frames for gene {gene_name}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        if seqs_with_frameshifts_dict[gene_name]:
            return target_file_basename, seqs_with_frameshifts_dict
        else:
            return target_file_basename


def check_and_correct_reading_frames_multiprocessing(new_target_sequences_folder, frameshifts_folder,
                                                     uncorrected_frameshifts_folder, exonerate_results_folder,
                                                     refs_for_exonerate, pool_threads,
                                                     skip_fix_frameshifts_with_exonerate=False):
    """
    Run check_and_correct_reading_frames for each new target sequence file via function
    <check_and_correct_reading_frames> via multiprocessing.
    """
    createfolder(frameshifts_folder)
    createfolder(uncorrected_frameshifts_folder)
    if not skip_fix_frameshifts_with_exonerate:
        createfolder(exonerate_results_folder)

    logger.info('===> Checking and correcting reading frames for sequences in new target files...')
    target_files = [file for file in sorted(glob.glob(f'{new_target_sequences_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(check_and_correct_reading_frames, target_file, frameshifts_folder,
                                      uncorrected_frameshifts_folder, exonerate_results_folder, refs_for_exonerate,
                                      counter, lock, num_files_to_process=len(target_files),
                                      skip_fix_frameshifts_with_exonerate=skip_fix_frameshifts_with_exonerate)
                          for target_file in target_files]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

        seqs_with_frameshifts = []
        for future in future_results:
            try:
                target_file, seqs_with_frameshifts_dict = future.result()
                seqs_with_frameshifts.append(seqs_with_frameshifts_dict)
            # except (ValueError, TypeError) as e:
            except ValueError:
                target_file = future.result()

    processed_target_files = [target_file for target_file in glob.glob(f'{new_target_sequences_folder}/*.fasta') if
                              file_exists_and_not_empty(target_file)]
    logger.info(f'\n{len(processed_target_files)} alignments processed from {len(future_results)} alignment files.\n')
    return seqs_with_frameshifts


def run_exonerate(seqs_with_frameshifts_dict, refs_for_exonerate, single_gene_new_target_file, gene_name,
                  gene_exonerate_folder):
    """
    For genes with stop codons in each forwards reading frame, writes files for Exonerate and performs and parses
    Exonerate output. Returns a dictionary of sequences that couldn't be fixed.

    Uses the function <correct_frameshifts>.
    """
    fixed_seqs = []
    try:
        commandline_refs_for_exonerate = list(sorted(set(refs_for_exonerate)))
    except TypeError:
        commandline_refs_for_exonerate = None
    if commandline_refs_for_exonerate:
        re_compile_string = '|'.join(commandline_refs_for_exonerate)
        logger.debug(f'Non-default references used for Exonerate searches...')
    else:
        re_compile_string = 'Ambtr|Arath|Orysa'
        logger.debug('Default references used for Exonerate searches')

    createfolder(gene_exonerate_folder)
    frameshift_free_seqs_list = list(SeqIO.parse(single_gene_new_target_file, 'fasta'))
    pattern = re.compile(re_compile_string)
    ref_names_seqs_and_lengths = [(seq.name, seq, len(seq.seq.ungap(gap='-'))) for seq in
                                  frameshift_free_seqs_list
                                  if re.match(pattern, seq.name)]
    longest_ref_name = max(ref_names_seqs_and_lengths, key=itemgetter(2))[0]
    longest_ref_seq_trans = pad_seq(max(ref_names_seqs_and_lengths, key=itemgetter(2))[1]).translate()
    prot_for_exonerate_name = f'{gene_exonerate_folder}/{gene_name}_{longest_ref_name}_protein.fasta'
    with open(prot_for_exonerate_name, 'w') as prot:
        SeqIO.write(longest_ref_seq_trans, prot, 'fasta')
    for seq in seqs_with_frameshifts_dict[gene_name]:
        frameshift_seq_for_exonerate_name = f'{gene_exonerate_folder}/{gene_name}_{seq.name}.fasta'
        exonerate_result_file = f'{gene_exonerate_folder}/{gene_name}_{seq.name}.exn'
        with open(frameshift_seq_for_exonerate_name, 'w') as frameshift_seq:
            SeqIO.write(seq, frameshift_seq, 'fasta')
        exonerate_command = f'exonerate -m protein2genome --revcomp False --showalignment yes ' \
                            f'--showvulgar no -V 0 --refine full {prot_for_exonerate_name} ' \
                            f'{frameshift_seq_for_exonerate_name} > {exonerate_result_file}'
        exonerate_command_norefine = f'exonerate -m protein2genome --revcomp False --showalignment yes ' \
                                     f'--showvulgar no -V 0 {prot_for_exonerate_name} ' \
                                     f'{frameshift_seq_for_exonerate_name} > {exonerate_result_file}'
        try:
            subprocess.run(exonerate_command, shell=True, check=True)
        except: # Sometimes using the '--refine full' option in Exonerate causes an error - presumably a bug?
            try:
                subprocess.run(exonerate_command_norefine, shell=True, check=True)
            except:
                logger.info(f'\nExonerate produced an error for gene {gene_name}, sequence {seq.name}. Command was: '
                            f'{exonerate_command_norefine}\n')
                continue
        if not file_exists_and_not_empty(exonerate_result_file):
            logger.info(f'\nNo Exonerate results for gene {gene_name}, sequence {seq.name}\n')
            continue
        fixed_seq = correct_frameshifts(exonerate_result_file)
        if fixed_seq:
            seqrecord = SeqRecord(Seq(fixed_seq), id=seq.id, name=seq.name,
                                  description='corrected_frameshifts')
            fixed_seqs.append(seqrecord)
        else:
            logger.debug(f"\nCouldn't fix frameshifts for gene {gene_name}, sequence {seq.name}, skipping\n")
            continue
    if len(fixed_seqs) != 0:
        fixed_seq_names = [seq.name for seq in fixed_seqs]
        unfixed_seqs = [seq for seq in seqs_with_frameshifts_dict[gene_name] if seq.name not in
                        fixed_seq_names]
        del seqs_with_frameshifts_dict[gene_name]
        for seq in unfixed_seqs:
            seqs_with_frameshifts_dict[gene_name].append(seq)
        fixed_seqs_file = f'{gene_exonerate_folder}/{gene_name}_fixed_seqs.fasta'
        with open(fixed_seqs_file, 'w') as fixed:
            SeqIO.write(fixed_seqs, fixed, 'fasta')
        with open(single_gene_new_target_file, 'w') as target_file_with_fixed_seqs:
            frameshift_free_seqs_list.extend(fixed_seqs)
            SeqIO.write(frameshift_free_seqs_list, target_file_with_fixed_seqs, 'fasta')
    return seqs_with_frameshifts_dict


def correct_frameshifts(exonerate_result_file):
    """
    Parses an exonerate results file using SearchIO and attempts to fix frameshifts. Checks if a full open reading
    frame can be found. If not, excludes the sequence.
    https://biopython.org/DIST/docs/api/Bio.SearchIO.ExonerateIO-module.html
    """

    all_qresult = list(SearchIO.parse(exonerate_result_file, 'exonerate-text'))
    for query_result in all_qresult:
        hit_name = query_result.items[0][0]
        for hit in query_result:  # Only a single hit in this case, with single query and single target?
            if len(hit.hsps) == 1:
                concatenated_single_hsp_alignment_seqs = []
                for alignment in hit.hsps[0].aln_annotation_all:
                    alignment_seq = ''.join(alignment['hit_annotation'])
                    concatenated_single_hsp_alignment_seqs.append(alignment_seq.replace('-', ''))
                query_range = hit.hsps[0].query_range_all
                hit_range = hit.hsps[0].hit_range_all

                # If there's only one hit range i.e no introns/frameshift:
                if len(hit_range) == 1:
                    concatenated_seq = str(concatenated_single_hsp_alignment_seqs[0])
                    num_stop_codons = pad_seq(Seq(concatenated_seq)).translate().count('*')
                    if num_stop_codons != 0:
                        logger.debug(f"Couldn't fix frameshifts for {hit_name}, skipping\n")
                        return
                    return concatenated_seq

                # If there more than one hit range (due to e.g. introns/frameshifts, insert Ns between seqs
                ns_to_insert_list = []
                for filtered_range_pair_hit in pairwise(query_range):  # Should this actually be hit range, instead? No, that would insert Ns where introns have been removed.
                    left_seq, right_seq = filtered_range_pair_hit
                    three_prime_upstream_seq = left_seq[-1]
                    five_prime_downstream_seq = right_seq[0]
                    codons_in_query_but_not_in_hit_seqs = five_prime_downstream_seq - three_prime_upstream_seq
                    num_ns_to_insert = codons_in_query_but_not_in_hit_seqs * 3
                    ns = 'N' * num_ns_to_insert
                    ns_to_insert_list.append(ns)
                zip_seqs_and_ns = [item for item in
                                   itertools.zip_longest(concatenated_single_hsp_alignment_seqs,
                                                         ns_to_insert_list, fillvalue='')]
                zip_seqs_and_ns_flattened = list(itertools.chain(*zip_seqs_and_ns))
                concatenated_seq = ''.join(zip_seqs_and_ns_flattened)
                num_stop_codons = pad_seq(Seq(concatenated_seq)).translate().count('*')
                if num_stop_codons != 0:
                    logger.debug(f"Couldn't fix frameshifts for {hit_name}, skipping\n")
                    return
                return concatenated_seq
            else: # i.e. If there's more than 1 hsp for the hit
                hsp_seq_dict = OrderedDict()
                query_ranges = []
                hit_ranges = []
                hsp_seq_number = 1
                for hsp in sorted(hit.hsps, key=lambda x: x.hit_end):
                    hsp_seq_name = f'contig_{hsp_seq_number}'
                    concatenated_hsp_alignment_seqs = []
                    for alignment in hsp.aln_annotation_all:
                        alignment_seq = ''.join(alignment['hit_annotation'])
                        concatenated_hsp_alignment_seqs.append(alignment_seq)
                    hsp_seq = ''.join(concatenated_hsp_alignment_seqs).replace('-', '')
                    hsp_seq_dict[hsp_seq_name] = hsp_seq
                    hit_ranges.append(hsp.hit_range_all)
                    query_ranges.append(hsp.query_range_all)
                    hsp_seq_number += 1
                hsp_seq_slices_dict = {key: [None, None] for key in hsp_seq_dict.keys()}
                filtered_ranges_hit = hit_ranges.copy()
                filtered_ranges_query = query_ranges.copy()
                contig_names = [key for key in hsp_seq_dict.keys()]

                range_comparisons = itertools.permutations(zip(hit_ranges, query_ranges, contig_names), 2)
                for combination in range_comparisons:
                    seq_a_hit_5prime_boundary = combination[0][0][0][0]
                    seq_a_hit_3prime_boundary = combination[0][0][-1][-1]
                    seq_a_hit_range = combination[0][0]
                    seq_a_query_range = combination[0][1]
                    seq_a_hit_name = combination[0][-1]
                    seq_b_hit_5prime_boundary = combination[1][0][0][0]
                    seq_b_hit_3prime_boundary = combination[1][0][-1][-1]
                    seq_b_hit_range = combination[1][0]

                    if seq_a_hit_5prime_boundary >= seq_b_hit_5prime_boundary and seq_a_hit_3prime_boundary \
                            <= seq_b_hit_3prime_boundary:
                        logger.debug(f'Range {seq_a_hit_range} is encapsulated is range {seq_b_hit_range}, '
                                     f'removing..,\n')
                        try:
                            filtered_ranges_hit.remove(seq_a_hit_range)
                            filtered_ranges_query.remove(seq_a_query_range)
                            del hsp_seq_dict[seq_a_hit_name]
                            del hsp_seq_slices_dict[seq_a_hit_name]
                        except ValueError:
                            logger.debug(f'Ranges already removed')
                            pass

                filtered_contig_name_list = [key for key in hsp_seq_dict.keys()]
                ns_to_insert_list = []
                for filtered_range_pair_hit, name_pair_hit, filtered_range_query_pair in \
                        zip(pairwise(filtered_ranges_hit), pairwise(filtered_contig_name_list),
                            pairwise(filtered_ranges_query)):

                    left_seq_hit_name, right_seq_hit_name = name_pair_hit
                    left_seq_hit_range, right_seq_hit_range = filtered_range_pair_hit
                    left_seq_query_range, right_seq_query_range = filtered_range_query_pair
                    three_prime_left_hit = left_seq_hit_range[-1][-1]
                    five_prime_right_hit = right_seq_hit_range[0][0]
                    three_prime_left_query = left_seq_query_range[-1][-1]
                    five_prime_right_query = right_seq_query_range[0][0]

                    if five_prime_right_hit < three_prime_left_hit:
                        ns_to_insert = ''
                        ns_to_insert_list.append(ns_to_insert)
                        overlap_offset_in_nucleotides = (three_prime_left_hit - five_prime_right_hit) - 1
                        five_prime_downstream_seq_slice = [overlap_offset_in_nucleotides, None]
                        hsp_seq_slices_dict[right_seq_hit_name] = five_prime_downstream_seq_slice
                    else:
                        codons_in_query_but_not_in_hit_seqs = five_prime_right_query - three_prime_left_query
                        num_ns_to_insert = codons_in_query_but_not_in_hit_seqs * 3
                        ns = 'N' * num_ns_to_insert
                        ns_to_insert_list.append(ns)
                        pass
                for key, value in hsp_seq_dict.items():
                    left_slice_coordinate = hsp_seq_slices_dict[key][0]
                    right_slice_coordinate = hsp_seq_slices_dict[key][1]
                    hsp_seq_dict[key] = hsp_seq_dict[key][left_slice_coordinate:right_slice_coordinate]
                zip_seqs_and_ns = [item for item in
                                   itertools.zip_longest(hsp_seq_dict.values(), ns_to_insert_list, fillvalue='')]
                zip_seqs_and_ns_flattened = list(itertools.chain(*zip_seqs_and_ns))
                concatenated_seq = ''.join(zip_seqs_and_ns_flattened)
                num_stop_codons = pad_seq(Seq(concatenated_seq)).translate().count('*')
                if num_stop_codons != 0:
                    logger.debug(f"Couldn't fix frameshifts for {hit_name}, skipping\n")
                    return
                return concatenated_seq


def grouped(iterable, n):
    """s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    Used in the function correct_frameshifts.
    """
    return zip(*[iter(iterable)]*n)


def pairwise(iterable):  # CJJ
    """s -> (s0,s1), (s1,s2), (s2, s3), ...
    Used in the function fullContigs to iterate over overlapping pairs of hit_start_and_end_indices.
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def create_mega_target_file(final_seqs_folder, outfolder):
    """
    Writes a new target file by concatenating the single gene fasta files.
    """
    createfolder(outfolder)
    megatarget_file = f'{outfolder}/BYO_target.fasta'
    if file_exists_and_not_empty(megatarget_file):
        os.remove(megatarget_file)
    final_seqs_for_megatarget = [file for file in glob.glob(f'{final_seqs_folder}/*.fasta')]
    concatenate_small(megatarget_file, *final_seqs_for_megatarget)
    logger.info(f"\n===> Concatenating {len(final_seqs_for_megatarget)} fasta files into the target file "
                f"'{os.path.basename(megatarget_file)}'")
    return megatarget_file


def megatarget_single_gene_alignments(final_seqs_file, output_folder, warnings_folder, collated_warnings, counter, lock,
                                      num_files_to_process, mafft_threads=2):
    """
    Create single genes alignments from fasta files.
    """
    basename = os.path.basename(final_seqs_file)
    expected_alignment_file = f'{output_folder}/{re.sub("[.]fasta", ".aln.fasta", basename)}'
    genes_with_warnings = [warning.split()[-1].strip('!') for warning in collated_warnings]

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f' Alignment exists for {basename}, skipping...')
    except AssertionError:
        # mafft_cline = (MafftCommandline(localpair=True, lep=-1.0, thread=mafft_threads, input=f'{final_seqs_file}'))
        mafft_cline = (MafftCommandline('linsi', thread=mafft_threads, input=f'{final_seqs_file}'))
        stdout, stderr = mafft_cline()
        with open(expected_alignment_file, 'w') as alignment_file:
            alignment_file.write(stdout)
        if basename.split('.')[0] in genes_with_warnings:
            shutil.copy(expected_alignment_file, warnings_folder)
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished aligning fasta file for {basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return expected_alignment_file


def megatarget_single_gene_alignments_multiprocessing(final_seqs_folder, output_folder, warnings_folder,
                                                      collated_warnings, mafft_threads=2, pool_threads=1):
    """
    Generate alignments via function <megatarget_single_gene_alignments> using multiprocessing.
    """
    createfolder(output_folder)
    createfolder(warnings_folder)
    logger.info('===> Creating alignments from final gene fasta files...')
    fasta_files = [file for file in sorted(glob.glob(f'{final_seqs_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(megatarget_single_gene_alignments, final_seqs_file, output_folder,
                                      warnings_folder, collated_warnings, counter, lock,
                                      num_files_to_process=len(fasta_files), mafft_threads=mafft_threads)
                          for final_seqs_file in fasta_files]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
        alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                          file_exists_and_not_empty(alignment)]
        logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} alignment files.')
    if len(fasta_files) != len(alignment_list):
        sys.exit(f'Only {len(alignment_list)} alignments were processed from {len(fasta_files)} genes, check '
                 f'for errors!')


def write_report(original_targetfile, transcriptome_folder, new_targetfile_folder, reports_folder,
                 long_seqs_warnings, seqs_with_frameshifts, seqs_with_ns):
    """
    Writes a report (printing to screen and logging to file) containing stats and warnings.
    """

    createfolder(reports_folder)

    transcriptome_names = [re.split('-|[.]', str(os.path.basename(name)))[0] for name in
                           glob.glob(f'{transcriptome_folder}/*renamed.fasta')]

    # Recover number of original sequences and gene names
    with open(original_targetfile) as original:
        gene_names = set()
        original_num_seqs = 0
        seqs = SeqIO.parse(original, 'fasta')
        for sequence in seqs:
            original_num_seqs += 1
            gene_names.add(sequence.name.split('-')[-1])
    gene_names_sorted = sorted(name for name in gene_names)

    # Recover number and names of transcriptome sequences added
    transcriptome_seqs_added = defaultdict(list)
    with open(f'{new_targetfile_folder}/BYO_target.fasta') as new:
        new_num_seqs = 0
        seqs = SeqIO.parse(new, 'fasta')
        for sequence in seqs:
            new_num_seqs += 1
            for name in transcriptome_names:
                if re.match(name, sequence.name):
                    transcriptome_seqs_added[name].append(sequence.name)

    # Get stats and write 'report_per_gene.csv' file
    summary_report = [['', 'Num_seqs_added', 'Num_seqs_grafted']]
    per_gene_report_final = []
    for transcriptome_name, gene_names in transcriptome_seqs_added.items():
        per_transcriptome_report = [transcriptome_name]
        for gene in gene_names_sorted:
            if gene in [gene_name.split('-')[-1] for gene_name in gene_names]:
                per_transcriptome_report.append('1')
            else:
                per_transcriptome_report.append('0')
        num_grafted = len([gene_name for gene_name in gene_names if re.search('grafted', gene_name)])
        per_gene_report_final.append(per_transcriptome_report)
        summary_report.append([transcriptome_name, str(len(gene_names)), str(num_grafted)])

    with open(f'{reports_folder}/report_per_gene.csv', 'w') as per_gene_outfile:
        per_gene_outfile.write(f",{','.join(gene_names_sorted)}\n")
        for line in per_gene_report_final:
            per_gene_outfile.write(f"{','.join(line)}\n")

    # Print summary report to screen and write to 'summary_report.csv' file
    logger.info(f'\n\n{"*" * 42} SUMMARY REPORT {"*" * 42}')
    logger.info(f"\n  Transcriptomes processed: {', '.join(transcriptome_names)}")
    logger.info(f'  Number of sequences in original target_file: {original_num_seqs}')
    logger.info(f'  Number of sequences in new target_file: {new_num_seqs}\n')
    with open(f'{reports_folder}/summary_report.csv', 'w') as summary:
        summary.write(f"Transcriptomes processed: {'; '.join(transcriptome_names)}\n")
        summary.write(f'Number of sequences in original target_file: {original_num_seqs}\n')
        summary.write(f'Number of sequences in new target_file: {new_num_seqs}\n\n')
        for entry in summary_report:
            logger.info(f'{entry[0]:>10}{entry[1]:>20}{entry[2]:>20}')
            summary.write(f"{','.join(entry)}\n")
        logger.info(f'\n')
        if long_seqs_warnings:
            logger.info(f'\n  ==== WARNING ====')

            fill = textwrap.fill(f'At least one newly added sequence is longer than expected for '
                                 f'{len(long_seqs_warnings)} genes. In most cases this is likely due to the presence '
                                 f'of an intron in the transcriptome sequence. Details have been logged to the summary '
                                 f'report, and the corresponding gene alignments have been written to folder '
                                 f'"16_final_seqs_alignments_with_warnings".', 98)
            logger.info(textwrap.indent(fill, '  '))
            summary.write(f'\n')
            for warning in long_seqs_warnings:
                summary.write(f'{warning}\n')
        else:
            logger.info(f'\n\n  No long sequences warnings!')
            summary.write(f'\nNo long sequence warnings!\n')
        if seqs_with_frameshifts:
            summary.write(f'\n')
            summary.write(f'GENE,SEQUENCES WITH FRAMESHIFTS (not included in new target file)\n')

            all_seqs_with_frameshift = []
            for gene_dict in seqs_with_frameshifts:
                for gene_name, seqs_with_frameshifts_list in gene_dict.items():
                    gene_names_iter = iter([seq.name for seq in seqs_with_frameshifts_list])
                    summary.write(f'{gene_name},{next(gene_names_iter)}\n')
                    for other_sequence in gene_names_iter:
                        summary.write(f',{other_sequence}\n')
                    simplified_name_fields = [re.split('_|-', seq.name) for seq in seqs_with_frameshifts_list]
                    seqs_with_frameshifts_names = [f'{gene_name_fields[-1]}-{gene_name_fields[0]}' for gene_name_fields
                                                   in simplified_name_fields]
                    all_seqs_with_frameshift.extend(seqs_with_frameshifts_names)
            logger.info(f'\n  ==== NOTE ====')
            fill = textwrap.fill(f'Some genes had sequences that could not be translated without stop codons; '
                                 f'{len(all_seqs_with_frameshift)} such sequences have been removed, affecting '
                                 f'{len(seqs_with_frameshifts)} genes. Details have been logged to the the summary '
                                 f'report, and the removed sequences have been written to folder '
                                 f'"13_gene_sequences_with_uncorrected_frameshifts".', 98)
            logger.info(textwrap.indent(fill, '  '))

        if seqs_with_ns:
            summary.write(f'\n')
            summary.write(f'GENE,SEQUENCES THAT CONTAIN(ED) Ns\n')
            all_seqs_with_ns = []
            for gene_dict in seqs_with_ns:
                for gene_name, seqs_with_ns_list in gene_dict.items():
                    seqs_with_ns_names = [f'{gene_name}-{seq.name.split("-")[1]}' for seq in seqs_with_ns_list]
                    seqs_with_ns_names_iter = iter(seqs_with_ns_names)
                    summary.write(f'{gene_name},{next(seqs_with_ns_names_iter)}\n')
                    for other_sequence in seqs_with_ns_names_iter:
                        summary.write(f',{other_sequence}\n')
                    all_seqs_with_ns.extend(seqs_with_ns_names)
            logger.info(f'\n  ==== NOTE ====')
            fill = textwrap.fill(f'Some genes had sequences that contain(ed) Ns; {len(all_seqs_with_ns)} such '
                                 f'sequences were found, affecting {len(seqs_with_ns)} genes. Details have been '
                                 f'logged to the the summary report, and the corresponding sequences have been written '
                                 f'to folder "07_gene_sequences_with_Ns".', 98)
            logger.info(textwrap.indent(fill, '  '))

        if seqs_with_ns and seqs_with_frameshifts:
            both_ns_and_frameshifts = list(set(all_seqs_with_ns) & set(all_seqs_with_frameshift))
            frameshift_but_no_ns = list(set(all_seqs_with_frameshift).difference(all_seqs_with_ns))

            if both_ns_and_frameshifts:
                logger.info(f'\n  ==== NOTE ====')
                both_ns_and_frameshifts_dict = defaultdict(list)
                for seq in both_ns_and_frameshifts:
                    gene_name = seq.split('-')[0]
                    both_ns_and_frameshifts_dict[gene_name].append(seq)
                fill = textwrap.fill(f'Some genes had sequences that contained both Ns and frameshifts that could not '
                                     f'be corrected; {len(both_ns_and_frameshifts)} such sequences were found. '
                                     f'Details have been logged to the the summary report.', 98)
                logger.info(textwrap.indent(fill, '  '))
                summary.write(f'\n')
                summary.write(f'GENE,SEQUENCES THAT CONTAINED BOTH FRAMESHIFTS AND Ns\n')
                for gene_name, seqs_list in both_ns_and_frameshifts_dict.items():
                    gene_names_iter = iter(seqs_list)
                    summary.write(f'{gene_name},{next(gene_names_iter)}\n')
                    for other_sequence in gene_names_iter:
                        summary.write(f',{other_sequence}\n')

            if frameshift_but_no_ns:
                frameshift_but_no_ns_dict = defaultdict(list)
                for seq in frameshift_but_no_ns:
                    gene_name = seq.split('-')[0]
                    frameshift_but_no_ns_dict[gene_name].append(seq)
                logger.info(f'\n  ==== NOTE ====')
                fill = textwrap.fill(f'Some genes had sequences that contained frameshifts but no Ns; '
                                     f'{len(frameshift_but_no_ns)} such sequences were found. Details have been logged '
                                     f'to the the summary report.', 98)
                logger.info(textwrap.indent(fill, '  '))
                summary.write(f'\n')
                summary.write(f'GENE,SEQUENCES THAT CONTAINED FRAMESHIFTS BUT NO Ns\n')
                for gene_name, seqs_list in frameshift_but_no_ns_dict.items():
                    seqs_list_iter = iter(seqs_list)
                    summary.write(f'{gene_name},{next(seqs_list_iter)}\n')
                    for other_sequence in seqs_list_iter:
                        summary.write(f',{other_sequence}')

    logger.info(f'\n{"*" * 100}\n')
    logger.info(textwrap.fill("This summary report has been written to '18_reports/summary_report.csv', and a per-gene "
                              "report has been written to '18_reports/report_per_gene.csv'",  98))
    return


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('target_file', type=str, help='target fasta file containing nucleotide sequences named '
                                                      'following the convention e.g. >AJFN-4471')
    parser.add_argument('transcriptomes_folder', type=str, help='Folder containing transcriptome files. These can '
                                                                'optionally be zipped in .zip, .gz or .bz2 format')
    parser.add_argument("-num_hits_to_recover", type=int, help='Number of HMM hits to recover from each transcriptome, '
                                                               'default is 1', default='1', metavar='<integer>')
    parser.add_argument("-python_threads", type=int, help='number of threads for multiprocessing pools, default is 1',
                        default='1', metavar='<integer>')
    parser.add_argument("-external_program_threads", type=int, help='number of threads for HMMer, mafft, default is '
                                                                    '1', default='1', metavar='<integer>')
    parser.add_argument("-length_percentage", type=float, help='length percentage cut-off for grafting or discarding, '
                                                               'default is 0.85', default='0.85', metavar='<float>')
    parser.add_argument("-hmmsearch_evalue", type=str, help='Evalue threshold for searching transcriptomes with HMM'
                                                            'profiles, default is 1e-50', default='1e-50',
                        metavar='<number in scientific notation; default is 1e-50>')
    parser.add_argument("-trim_to_refs", action='append', type=str, dest='refs_for_manual_trimming',
                        help='Name corresponding to a taxon present in target file, used for trimming and correcting '
                             'frameshifts in transcriptome hits', metavar='<taxon_name>')
    parser.add_argument("-no_n", dest="no_n", action='store_true', help="Remove n characters from transcriptomes hits",
                        default=False)
    parser.add_argument("-skip_exonerate_frameshift_fix", dest="no_exonerate_fix", action='store_true',
                        help="Do not use Exonerate to fix frameshifts; discard such sequences instead", default=False)
    # parser.print_help()
    results = parser.parse_args()
    return results


########################################################################################################################
########################################################################################################################

def main():
    results = parse_arguments()

    folder_01 = cwd + '/01_target_gene_fasta_files'
    folder_02 = cwd + '/02_target_gene_alignment'
    folder_03 = cwd + '/03_target_gene_hmms'
    folder_04 = cwd + '/04_target_gene_hmm_hits'
    folder_05 = cwd + '/05_target_gene_hmm_logs'
    folder_06 = cwd + '/06_align_extractions'
    folder_07 = cwd + '/07_gene_sequences_with_Ns'
    folder_08 = cwd + '/08_trimmed_alignments'
    folder_09 = cwd + '/09_grafting_sequence_alignments'
    folder_10 = cwd + '/10_seqs_to_add'
    folder_11 = cwd + '/11_final_seqs'
    folder_12 = cwd + '/12_gene_sequences_with_frameshifts'
    folder_13 = cwd + '/13_gene_sequences_with_uncorrected_frameshifts'
    folder_14 = cwd + '/14_gene_exonerate_results'
    folder_15 = cwd + '/15_final_seqs_alignments'
    folder_16 = cwd + '/16_final_seqs_alignments_with_warnings'
    folder_17 = cwd + '/17_mega_target_file'
    folder_18 = cwd + '/18_reports'

    check_dependencies(results.target_file,
                       results.transcriptomes_folder,
                       results.python_threads,
                       results.external_program_threads,
                       results.length_percentage,
                       results.hmmsearch_evalue,
                       results.refs_for_manual_trimming,
                       results.no_n)

    check_files_for_processing(results.target_file,
                               results.transcriptomes_folder,
                               results.refs_for_manual_trimming)

    split_targets(results.target_file,
                  output_folder=folder_01)

    align_targets_multiprocessing(target_gene_folder=folder_01,
                                  alignments_output_folder=folder_02,
                                  algorithm='linsi',
                                  pool_threads=results.python_threads,
                                  mafft_threads=results.external_program_threads)

    create_hmm_profile_multiprocessing(alignments_folder=folder_02,
                                       hmm_output_folder=folder_03,
                                       pool_threads=results.python_threads,
                                       hmmbuild_threads=results.external_program_threads)

    hmm_vs_transcriptome_multiprocessing(hmmprofile_folder=folder_03,
                                         transcriptomes_folder=results.transcriptomes_folder,
                                         hits_output_folder=folder_04, hmm_logs_output_folder=folder_05,
                                         num_hits_to_recover=results.num_hits_to_recover,
                                         pool_threads=results.python_threads,
                                         hmmsearch_threads=results.external_program_threads,
                                         hmmsearch_evalue=results.hmmsearch_evalue)

    seqs_with_ns = align_extractions_multiprocessing(
        alignments_folder=folder_02,
        output_folder=folder_06,
        hit_folder=folder_04,
        seqs_with_ns_folder=folder_07,
        concatenated_folder='concatenated_hits.tmp',
        mafft_threads=results.external_program_threads,
        pool_threads=results.python_threads,
        no_n=results.no_n)


    long_seqs_warnings = trim_and_discard_or_graft_multiprocessing(alignments_folder=folder_06,
                                                                   trimmed_alignment_folder=folder_08,
                                                                   alignments_for_grafting_folder=folder_09,
                                                                   new_sequence_folder=folder_10,
                                                                   pool_threads=results.python_threads,
                                                                   mafft_threads=results.external_program_threads,
                                                                   length_percentage=results.length_percentage,
                                                                   refs_for_trimmed=results.refs_for_manual_trimming,
                                                                   discard_short=False)

    fasta_files = [file for file in glob.glob(f'{folder_01}/*.fasta')]
    for counter, fasta_file in enumerate(fasta_files, 1):
        create_new_targets_file(fasta_file,
                                seqs_to_add_folder=folder_10,
                                seqs_to_write_folder=folder_11)

    seqs_with_frameshifts = check_and_correct_reading_frames_multiprocessing(
        new_target_sequences_folder=folder_11,
        frameshifts_folder=folder_12,
        uncorrected_frameshifts_folder=folder_13,
        exonerate_results_folder=folder_14,
        refs_for_exonerate=results.refs_for_manual_trimming,
        pool_threads=results.python_threads,
        skip_fix_frameshifts_with_exonerate=results.no_exonerate_fix)


    megatarget_single_gene_alignments_multiprocessing(final_seqs_folder=folder_11,
                                                      output_folder=folder_15,
                                                      warnings_folder=folder_16,
                                                      collated_warnings=long_seqs_warnings,
                                                      mafft_threads=results.external_program_threads,
                                                      pool_threads=results.python_threads)

    create_mega_target_file(final_seqs_folder=folder_11,
                            outfolder=folder_17)

    write_report(results.target_file,
                 transcriptome_folder=results.transcriptomes_folder,
                 new_targetfile_folder=folder_17,
                 reports_folder=folder_18,
                 long_seqs_warnings=long_seqs_warnings,
                 seqs_with_frameshifts=seqs_with_frameshifts,
                 seqs_with_ns=seqs_with_ns)


########################################################################################################################
########################################################################################################################

if __name__ == '__main__':
    sys.exit(main())

########################################################################################################################
########################################################################################################################
