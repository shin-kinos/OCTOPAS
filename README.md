# ThesisProject

## Description 

* A Novel Algorithm to Reduce Redundancy of Large Size Phylogenetic Datasets

## Dependencies

**Test environment:**

* macOS Monterey (v12.6)
* Crescent server
* Python (v3.9.9)
* phylotreelib (v1.18.1)

## Installation

This script requires `phylotreelib` package (https://pypi.org/project/phylotreelib/).

So, if necessary, install this package by using `pip` command:

```
% pip install phylotreelib
```

Then, download the script `script_prototype.py` in `script` directory in your environment.

You can also download a small test phylogenetic tree data `test_data_500leaves.newick` in `test_data` direcotry, which is composed of 500 leaves of Human hemoglobin subunit alpha AA sequence and its homologous sequences.

## Usage

[e.g.]

```
% python3 script_prototype.py -i test_data_500leaves.newick -o output -O outdir
```

type `-h` or `--help` to see all the options supposed:

```
Options :

-i | --input   <String>     Input file name, REQUIRED
-o | --output  <String>     Output file name prefix, REQUIRED
-O | --outdir  <String>     Output directory name, if it is empty, all output files are saved in current directory
-f | --format  <String>     Format of input file:
                                * newick = Newick (default)
                                * nexus  = NEXUS
-n | --NL      <Integer>    Number of leaves remain, more than 0, default 3
-r | --RTL     <Real>       Threshold of relative tree length, at range of [0, 1], default 0.95
-s | --silent  <Flag>       If true, program does not show any logs on terminal during calculation
-v | --version <Flag>       Print current version, ignore all other arguments
-h | --help    <Flag>       Print this help, ignore all other arguments
```
