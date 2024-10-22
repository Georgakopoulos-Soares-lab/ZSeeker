ZDNA tool
==============

# Installation
```bash
pip install zdna_backend
```

# CLI Usage
```bash
zdna_backend --path ./test_GCA_f.fasta --n_jobs 1 --method=coverage
```

# Example: In Code usage
```python
from zdna_calculator import ZDNACalculatorSeq, Params
# Define parameters
params = Params(
    GC_weight=1.0,
    AT_weight=0.5,
    GT_weight=0.3,
    AC_weight=0.2,
    mismatch_penalty_starting_value=5,
    mismatch_penalty_linear_delta=2,
    mismatch_penalty_type='linear',
    method='coverage',
    threshold=10,
    consecutive_AT_scoring=[1, 2, 2],
    display_sequence_score=1
)

# Create a ZDNACalculatorSeq instance
zdna_calculator = ZDNACalculatorSeq(data="ACGTACGTACGT", params=params)

# Calculate subarrays above threshold
subarrays = zdna_calculator.subarrays_above_threshold()

# Print results
print(subarrays)
```

# Command-line Help
```bash
usage: zdna_backend.py [-h] [--path PATH] [--GC_weight GC_WEIGHT]
                       [--AT_weight AT_WEIGHT] [--GT_weight GT_WEIGHT]
                       [--AC_weight AC_WEIGHT]
                       [--mismatch_penalty_starting_value MISMATCH_PENALTY_STARTING_VALUE]
                       [--mismatch_penalty_linear_delta MISMATCH_PENALTY_LINEAR_DELTA]
                       [--mismatch_penalty_type {linear,exponential}]
                       [--method {coverage,score}] [--n_jobs N_JOBS]
                       [--threshold THRESHOLD]
                       [--consecutive_AT_scoring CONSECUTIVE_AT_SCORING]
                       [--max_resources_threshold MAX_RESOURCES_THRESHOLD]
                       [--display_sequence_score {0,1}]

Given a fasta file and the corresponding parameters it calculates the ZDNA
for each sequence present.

optional arguments:
  -h, --help            show this help message and exit
  --path PATH           Path to file analyzed
  --GC_weight GC_WEIGHT
                        Weight given to GC and CG transitions.
  --AT_weight AT_WEIGHT
                        Weight given to AT and TA transitions.
  --GT_weight GT_WEIGHT
                        Weight given to GT and TG transitions.
  --AC_weight AC_WEIGHT
                        Weight given to AC and CA transitions.
  --mismatch_penalty_starting_value MISMATCH_PENALTY_STARTING_VALUE
                        Penalty applied to the first non purine/pyrimidine
                        transition encountered.
  --mismatch_penalty_linear_delta MISMATCH_PENALTY_LINEAR_DELTA
                        Determines the rate of increase of the penalty for
                        every subsequent non purine/pyrimidine transition.
  --mismatch_penalty_type {linear,exponential}
                        Method of scaling the penalty for contiguous non
                        purine/pyrimidine transitions.
  --method {coverage,score}
                        Method used for the Z-DNA scoring algorithm.
  --n_jobs N_JOBS       Number of threads to use. Defaults to -1, which uses
                        the maximum available threads on CPU.
  --threshold THRESHOLD
                        Scoring threshold for a sequence to be considered
                        potentially Z-DNA forming.
  --consecutive_AT_scoring CONSECUTIVE_AT_SCORING
                        Penalty array for consecutive AT repeats forming
                        hairpin structures.
  --max_resources_threshold MAX_RESOURCES_THRESHOLD
                        Maximum resources threshold.
  --display_sequence_score {0,1}
                        Display the total sequence score (1) or not (0).
```
