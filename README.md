# ZSeeker

==============

# Installation
```bash
pip install ZSeeker
```

# CLI Usage
```bash
ZSeeker --fasta ./test_GCA_f.fasta --n_jobs 1
```

# Example: In Code usage
```python
from zseeker.zdna_calculator import ZDNACalculatorSeq, Params
# Define parameters
params = Params(
    GC_weight=5.0,
    AT_weight=0.5,
    GT_weight=1.1,
    AC_weight=1.3,
    mismatch_penalty_starting_value=4,
    mismatch_penalty_linear_delta=2,
    mismatch_penalty_type='linear',
    threshold=10,
    consecutive_AT_scoring=[1, 2, 2],
    display_sequence_score=1
    drop_threshold=50,
    total_sequence_scoring=False
)

# Create a ZDNACalculatorSeq instance and nput sequence
zdna_calculator = ZDNACalculatorSeq(data="ACGTACGTACGT", params=params)

# Calculate subarrays above threshold
subarrays = zdna_calculator.subarrays_above_threshold()

# Print results
print(subarrays)
```

# Command-line Help
```bash
usage: ZSeeker [-h] [--fasta FASTA] [--GC_weight GC_WEIGHT]
               [--AT_weight AT_WEIGHT] [--GT_weight GT_WEIGHT]
               [--AC_weight AC_WEIGHT]
               [--mismatch_penalty_starting_value MISMATCH_PENALTY_STARTING_VALUE]
               [--mismatch_penalty_linear_delta MISMATCH_PENALTY_LINEAR_DELTA]
               [--mismatch_penalty_type {linear,exponential}]
               [--n_jobs N_JOBS] [--threshold THRESHOLD]
               [--consecutive_AT_scoring CONSECUTIVE_AT_SCORING]
               [--display_sequence_score {0,1}]
               [--output_dir OUTPUT_DIR]
               [--gff_file GFF_FILE]
               [--drop_threshold DROP_THRESHOLD]
               [--total_sequence_scoring]

Given a fasta file and the corresponding parameters it calculates the
ZDNA for each sequence present.

options:
  -h, --help            show this help message and exit
  --fasta FASTA         Path to file analyzed
  --GC_weight GC_WEIGHT
                        Weight given to GC and CG transitions.
                        Default = 7.0
  --AT_weight AT_WEIGHT
                        Weight given to AT and TA transitions.
                        Default = 0.5
  --GT_weight GT_WEIGHT
                        Weight given to GT and TG transitions.
                        Default = 1.25
  --AC_weight AC_WEIGHT
                        Weight given to AC and CA transitions.
                        Default = 1.25
  --mismatch_penalty_starting_value MISMATCH_PENALTY_STARTING_VALUE
                        Penalty applied to the first non
                        purine/pyrimidine transition encountered.
                        Default = 3
  --mismatch_penalty_linear_delta MISMATCH_PENALTY_LINEAR_DELTA
                        Only applies if penalty type is set to
                        linear. Determines the rate of increase of
                        the penalty for every subsequent non
                        purine/pyrimidine transition. Default = 3
  --mismatch_penalty_type {linear,exponential}
                        Method of scaling the penalty for contiguous
                        non purine/pyrimidine transition. Default =
                        linear
  --n_jobs N_JOBS       Number of threads to use. Defaults to -1,
                        which uses the maximum available threads on
                        CPU
  --threshold THRESHOLD
                        Scoring threshold for a for a sequence to be
                        considered potentially Z-DNA forming and
                        returned by the program. This parameter is
                        also used for determining how big the scoring
                        drop within a sequence should be, before it
                        is split into two separate Z-DNA candidate
                        sequences. Default=50
  --consecutive_AT_scoring CONSECUTIVE_AT_SCORING
                        Consecutive AT repeats form a hairpin
                        structure instead of Z-DNA. In order to
                        reflect that, a penalty array is defined,
                        which provides the score adjustment for the
                        first and the subsequent TA appearances. The
                        last element will be applied to every
                        subsequent TA appearance. For more
                        information see documentation. Default =
                        (0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0)
  --display_sequence_score {0,1}
  --output_dir OUTPUT_DIR
  --gff_file GFF_FILE Optional GFF file for gene annotation. Only 'gene' features are used.
  --drop_threshold DROP_THRESHOLD
                        Drop threshold used within subarrays detection logic. Acts as earlier stopping threshold. Lower values result in smaller Z-DNA sequences and larger values result in fewer but larger Z-DNA sequences.
  --total_sequence_scoring
                        If set, calculate the total score of all provided sequences, without looking for Z-DNA subsequences. Useful for researchers who have short sequences and want to estimate their Z-DNA potential.
```


# Example output file

|Chromosome|Start|End  |Z-DNA Score|Sequence                                    |
|----------|-----|-----|-----------|--------------------------------------------|
|Z1        |0.0  |15.0 |87.0       |TGCGTGCGCGCGCGCG                            |
|Z2        |0.0  |15.0 |87.0       |GCGCCCGCGCGCGCGC                            |
|Z3        |0.0  |11.0 |71.0       |GCGCGCGCGCGT                                |
|Z4        |0.0  |11.0 |65.0       |GCGCGTGCGCGC                                |
|Z5        |0.0  |10.0 |70.0       |CGCGCGCGCGC                                 |
|Z6        |0.0  |15.0 |63.0       |GCACGCACACGCGCGT                            |
|Z7        |0.0  |10.0 |70.0       |GCGCGCGCGCG                                 |
|Z8        |0.0  |13.0 |61.0       |CGCACGCGCACGCA                              |
|Z9        |0.0  |11.0 |59.0       |CGCGCGCGCACA                                |


# Example output file with annotations

|Chromosome|Start|End  |Z-DNA Score|Sequence                                    |gene_start|gene_end|gene_id       |gene_biotype  |strand|distance|distance_from_TSS|distance_from_TES|
|----------|-----|-----|-----------|--------------------------------------------|----------|--------|--------------|--------------|------|--------|-----------------|-----------------|
|AE004438.1|364  |391  |63.0       |ACGGTGCCGCAGCGGCCGTGTCGCCAGC                |362       |812     |gene-VNG_6001H|protein_coding|-     |0       |420              |2                |
|AE004438.1|2317 |2335 |51.5       |GCGGCGAGTCGCCGTCGCG                         |1904      |3719    |gene-VNG_6007H|protein_coding|-     |0       |1383             |413              |
|AE004438.1|3528 |3538 |52.75      |ACGTGCGCGCG                                 |1904      |3719    |gene-VNG_6007H|protein_coding|-     |0       |180              |1624             |
|AE004438.1|12771|12814|109.25     |GCTGTCGCTGTCGGCGGCGGCTGCCGCCGACGCGACAGCGTCGC|12846     |13380   |gene-VNG_6015H|protein_coding|-     |32      |565              |32               |
|AE004438.1|13178|13195|56.0       |ACGGCGCGTCAGCGGCGT                          |12846     |13380   |gene-VNG_6015H|protein_coding|-     |0       |184              |332              |
|AE004438.1|13533|13552|52.75      |ACGGCGCACCGCCAGCGTGT                        |12846     |13380   |gene-VNG_6015H|protein_coding|-     |153     |154              |687              |
|AE004438.1|13853|13872|70.0       |CGTCGGCGCACGCGCCGACG                        |14307     |15582   |gene-VNG_6016H|protein_coding|+     |435     |435              |1709             |
|AE004438.1|14960|14971|51.25      |GCGCGGTCGCGC                                |14307     |15582   |gene-VNG_6016H|protein_coding|+     |0       |653              |610              |
|AE004438.1|15105|15126|61.0       |CGCGTCGTCGGCGTCCGCGACG                      |14307     |15582   |gene-VNG_6016H|protein_coding|+     |0       |798              |455              |

# ZSeeker web application

The web version of ZSeeker can be found at:
[ZSeeker web application](https://zseeker.netlify.app)

And a dockerized version of it can be found at this repository for local deployments:
[ZSeeker web application dockerized](https://github.com/Georgakopoulos-Soares-lab/ZSeeker_docker)

