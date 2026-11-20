# 🧬 DNA Sequence Analyzer

A command-line bioinformatics toolkit written in **pure Python**, with optional **Biopython** integration to cross-check results against a trusted reference library.

Built as a learning project to implement core sequence-analysis algorithms from scratch.

## Features

- **FASTA parsing** (multi-record) and sequence validation
- **Composition**: length, base counts, GC content, GC content in sliding windows
- **Complement, reverse complement, transcription (DNA → RNA)**
- **Translation** to protein in any of 3 frames (standard genetic code)
- **Motif search** with IUPAC ambiguity codes, plus common restriction sites (EcoRI, BamHI, HindIII, NotI, PstI)
- **ORF finder** across all 6 reading frames
- **Pairwise alignment**: Needleman-Wunsch (global) and Smith-Waterman (local)
- **Charts**: base composition and GC-content plots (matplotlib)
- **Biopython cross-check**: verifies GC content, reverse complement, translation and alignment scores
- Unit tests with `pytest`

## Installation

```bash
git clone https://github.com/<your-username>/dna-sequence-analyzer.git
cd dna-sequence-analyzer
pip install -r requirements.txt   # optional: only needed for plots, Biopython and tests
```

The core toolkit needs **no dependencies**, only Python 3.8+.

## Usage

```bash
# Full summary of every sequence in a FASTA file
python main.py analyze --file examples/sample.fasta

# Or pass a sequence directly
python main.py analyze --seq ATGGCCATTGTAATGGGCCGCTGA

# Translate (frame 1-3), stopping at the first stop codon
python main.py translate --seq ATGGCCATTGTAATGGGCCGCTGA --frame 1 --to-stop

# Find motifs / restriction sites
python main.py motif --file examples/sample.fasta --restriction
python main.py motif --seq GACTCGATTC --motif GANTC

# Find ORFs of at least 30 amino acids
python main.py orfs --file examples/sample.fasta --min-aa 30

# Align two sequences
python main.py align GATTACA GCATGCT
python main.py align TTTTGATTACATTTT GATTACA --mode local
python main.py align GATTACA GCATGCT --compare-biopython

# Charts and Biopython cross-check
python main.py plot --file examples/sample.fasta --out summary.png
python main.py check --file examples/sample.fasta
```

### Example output

```
$ python main.py analyze --file examples/sample.fasta

== sample_gene_2 synthetic example ==
Length            : 76 bp
Base counts       : {'A': 19, 'C': 18, 'G': 22, 'T': 17}
GC content        : 52.63%
Protein (frame 1) : MKFGPKFGP*GSGLQITASARYMRG

$ python main.py align GATTACA GCATGCT
GATTACA
|  | | 
GCATGCT
Score: -1 | Identity: 42.9%
```

## Use as a library

```python
from dna_analyzer import gc_content, reverse_complement, translate, find_orfs, needleman_wunsch

seq = "ATGGCCATTGTAATGGGCCGCTGA"
print(gc_content(seq))
print(translate(seq, to_stop=True))
print(needleman_wunsch("GATTACA", "GCATGCT"))
```

## Project structure

```
dna_analyzer/
├── sequence.py         # FASTA reading, validation, GC content, complement
├── translate.py        # transcription and translation (codon table)
├── motifs.py           # motif search and ORF finder
├── align.py            # Needleman-Wunsch and Smith-Waterman
├── biopython_tools.py  # optional Biopython cross-checks
├── plot.py             # matplotlib charts
└── cli.py              # command-line interface
tests/                  # pytest unit tests
examples/sample.fasta   # synthetic example data
```

## Running tests

```bash
pytest
```

## How the algorithms work

- **GC content**: percentage of G and C bases; ambiguous `N` bases are ignored.
- **ORF finding**: scans each of the 6 frames for `ATG` and reads codons until a stop codon (`TAA`, `TAG`, `TGA`).
- **Needleman-Wunsch**: dynamic programming over a score matrix, then traceback from the bottom-right cell for the best end-to-end alignment.
- **Smith-Waterman**: same idea, but scores never drop below 0 and the traceback starts at the highest score, giving the best *local* match.

Default scoring: match `+1`, mismatch `-1`, gap `-2` (adjustable with `--match`, `--mismatch`, `--gap`).

## Ideas for future work

- Affine gap penalties and substitution matrices
- Streamlit web interface
- Codon usage tables and GC-skew analysis
- Phylogenetic tree construction from multiple sequences

## License

MIT. See [LICENSE](LICENSE).
