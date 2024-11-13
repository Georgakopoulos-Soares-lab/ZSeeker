from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
from pathlib import Path
import pandas as pd
import concurrent.futures
import multiprocessing as mp
import os
from termcolor import colored
import time
from typing import Callable, Iterator
import numpy as np
import functools
from .zdna_calculator import ZDNACalculatorSeq, Params
from collections.abc import Iterable
import logging
from attrs import define, field
import argparse

@define(kw_only=True, slots=True, frozen=True)
class ZDNASubmissionForm:
    recordSeq: str = field(converter=str)
    recordID: str = field(converter=str)
    input_fasta: os.PathLike[str] = field(converter=str)

def timeit(func: Callable) -> Callable:
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(colored("ZDNA extraction initialized.", "magenta"))
        then = time.perf_counter()
        result = func(*args, **kwargs)
        now = time.perf_counter()
        print(colored(f"Process finished within {now-then:.2f} second(s).", "magenta"))
        return result
    return wrapper

def read_fasta(fasta: Path | str) -> Iterator[tuple[str, str]]:
    if Path(fasta).name.endswith(".gz"):
        file = gzip.open(fasta, 'rt')
    else:
        file = open(fasta, encoding='utf-8', mode='r')
    for record in SimpleFastaParser(file):
        yield record[0].split(' ')[0], record[1]
    file.close()

def _extract(ID: int, sequence: str, params: Params) -> pd.DataFrame:
    zdna_string = ZDNACalculatorSeq(data=sequence, params=params)
    subarrays_detected = zdna_string.subarrays_above_threshold()
    subarrays_detected = pd.DataFrame(subarrays_detected, columns=params.headers[1:-1])
    scoring_array = zdna_string.scoring_array
    
    # perhaps delete?
    if not len(subarrays_detected):
        null_subarrays = pd.DataFrame([[np.nan]*len(params.headers)], columns=params.headers)
        # Explicitly set the 'Chromosome' column to object dtype to avoid FutureWarning
        null_subarrays['Chromosome'] = null_subarrays['Chromosome'].astype(object)
        null_subarrays.loc[:, "Chromosome"] = ID
        null_subarrays.loc[:, "totalSequenceScore"] = sum(scoring_array)
        return null_subarrays
    # <<<<<<<<<<<<<<<<<
    subarrays_detected.loc[:, "Chromosome"] = ID
    subarrays_detected.loc[:, "totalSequenceScore"] = sum(scoring_array)
    return subarrays_detected[params.headers]

def _subextraction(submission_forms: list[ZDNASubmissionForm], params: Params) -> pd.DataFrame:
    # process concurrently
    outputs = []
    for form in submission_forms:
        try:
            present_df = _extract(ID=form.recordID,
                                sequence=form.recordSeq,
                                params=params
                                )
            outputs.append(present_df)
        except IndexError as e:
            logging.error(
                f"Failed to process input fasta '{Path(form.input_fasta).name}' on chromosome {str(form.recordID)}. "
                f"Chromosome length: {len(str(form.recordSeq))}."
            )
            logging.error(f"Error message: '{e}'.")
            continue

    if outputs:
        zdna_df = pd.concat(outputs, axis=0).reset_index(drop=True)
    else:
        zdna_df = pd.DataFrame(columns=params.headers)
        logging.info("Empty dataframe derived from subprocess _subextraction")
    return zdna_df

def assign_tasks(tasks: list, total_buckets: int) -> list[list]:
    then = time.perf_counter()
    total = len(tasks)
    step = total // total_buckets
    remainder = total % total_buckets
    assigned_tasks = []
    infimum = 0
    while True:
        if remainder > 0:
            supremum = infimum + step + 1
            remainder -= 1
        else:
            supremum = infimum + step
        assigned_tasks.append(tasks[infimum: supremum])
        if len(assigned_tasks) == total_buckets:
            break
        infimum = supremum
    now = time.perf_counter()
    print(colored(f"Task assignment completed within {now-then:.2f} second(s).", "green"))
    return assigned_tasks

def extract_zdna_v2(fasta: str | os.PathLike[str], params: Params) -> pd.DataFrame:
    then = time.perf_counter()
    zdna_df = pd.DataFrame(columns=params.headers)
    outputs: list[pd.DataFrame] = []
    n_jobs = params.n_jobs
    assert isinstance(n_jobs, int) and n_jobs > 0, "Number of jobs must be a positive int."
    submission_forms: list[ZDNASubmissionForm] = []
    for rec_id, rec_seq in read_fasta(fasta):
        submission_forms.append(
            ZDNASubmissionForm(
                recordSeq=rec_seq.upper(),
                recordID=rec_id,
                input_fasta=fasta
            )
        )
    total_submission_forms = len(submission_forms)
    
    # split into records
    assigned_tasks = assign_tasks(submission_forms, n_jobs)
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_jobs,
                                                mp_context=mp.get_context("spawn")
                                                ) as executor:
        results = executor.map(_subextraction, assigned_tasks, [params]*n_jobs)
        for result in results:
            print("Task finished.")
            outputs.append(result)

    zdna_df = pd.concat(outputs, axis=0).reset_index(drop=True)
    seqID_observed = zdna_df['Chromosome'].nunique()
    seqID_without_zdna = zdna_df[zdna_df['Start'].isna()]['Chromosome']
    if seqID_without_zdna.shape[0] > 0:
        logging.warning(f"The following {seqID_without_zdna.shape[0]} sequence IDs were found without Z-DNA")
        for seqID in seqID_without_zdna:
            print(seqID)

    if total_submission_forms > seqID_observed:
        logging.warning("Z-DNA was not detected in all submitted sequence IDs.")

    now = time.perf_counter()
    print(f"Process finished within {now-then:.2f} second(s).")
    if params.display_sequence_score == 0:
        zdna_df.drop(columns=['totalSequenceScore'], inplace=True)
        zdna_df.dropna(inplace=True, axis=0)
    return zdna_df

def parse_consecutive_AT_scoring(value: str) -> tuple[float, ...]:
    # Split the string by commas and convert each value to float
    try:
        # Strip whitespace and split by comma
        values = [float(x.strip()) for x in value.split(',')]
        return tuple(values)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Invalid consecutive_AT_scoring format: {e}")


@timeit
def transform(path: Path | str, params: Params) -> pd.DataFrame:
    logging.info(f"Processing file '{Path(path).name}'...")
    output_dir = params.output_dir
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    # lazy loading
    zdna_df = extract_zdna_v2(path, params)
    print(colored("Saving Z-DNA dataframe...", "magenta"))
    output_file = output_dir.joinpath(f"{Path(path).stem}_zdna_score.csv")
    zdna_df.to_csv(output_file, mode="w", index=False)
    logging.info(f"File '{Path(path).name}' has been processed successfully.")
    return zdna_df

def main():
    parser = argparse.ArgumentParser(description="""Given a fasta file and the corresponding parameters
                                                    it calculates the ZDNA for each
                                                    sequence present.""")
    parser.add_argument("--fasta", type=str, default="test_file.fna", help="Path to file analyzed")
    parser.add_argument("--GC_weight", type=float, default=Params.GC_weight,
                        help=f"Weight given to GC and CG transitions. Default = {Params.GC_weight}")
    parser.add_argument("--AT_weight", type=float, default=Params.AT_weight,
                        help=f"Weight given to AT and TA transitions. Default = {Params.AT_weight}")
    parser.add_argument("--GT_weight", type=float, default=Params.GT_weight,
                        help=f"Weight given to GT and TG transitions. Default = {Params.GT_weight}")
    parser.add_argument("--AC_weight", type=float, default=Params.AC_weight,
                        help=f"Weight given to AC and CA transitions. Default = {Params.AC_weight}")
    parser.add_argument("--mismatch_penalty_starting_value", type=int, default=Params.mismatch_penalty_starting_value,
                        help=f"Penalty applied to the first non purine/pyrimidine transition encountered. "
                             f"Default = {Params.mismatch_penalty_starting_value}")
    parser.add_argument("--mismatch_penalty_linear_delta", type=int, default=Params.mismatch_penalty_linear_delta,
                        help=f"Only applies if penalty type is set to linear. Determines the rate of increase of the "
                             f"penalty for every subsequent non purine/pyrimidine transition. "
                             f"Default = {Params.mismatch_penalty_linear_delta}")
    parser.add_argument("--mismatch_penalty_type", choices=Params.mismatch_penalty_choices,
                        default=Params.mismatch_penalty_type,
                        help=f"Method of scaling the penalty for contiguous non purine/pyrimidine transition. "
                             f"Default = {Params.mismatch_penalty_type}")
    parser.add_argument("--n_jobs", type=int, default=Params.n_jobs,
                        help="Number of threads to use. Defaults to -1, which uses the maximum available threads on CPU")
    parser.add_argument("--threshold", type=int, default=Params.threshold,
                        help=f"Scoring threshold for a for a sequence to be considered potentially Z-DNA forming and "
                             f"returned by the program. This parameter is also used for determining how big the scoring "
                             f"drop within a sequence should be, before it is split into two separate Z-DNA candidate "
                             f"sequences. Default={Params.threshold}")
    parser.add_argument( "--consecutive_AT_scoring", type=parse_consecutive_AT_scoring, default=Params.consecutive_AT_scoring,
        help=f"Consecutive AT repeats form a hairpin structure instead of Z-DNA. In order to reflect that, "
            f"a penalty array is defined, which provides the score adjustment for the first and the "
            f"subsequent TA appearances. The last element will be applied to every subsequent TA "
            f"appearance. For more information see documentation. Default = {Params.consecutive_AT_scoring}"
    )
    parser.add_argument("--display_sequence_score", type=int, choices=[0, 1], default=0)
    parser.add_argument("--output_dir", type=str, default="zdna_extractions")

    args = parser.parse_args()
    _params = vars(args)
    fasta = _params.pop("fasta")
    # validation
    for key in _params:
        if key in ("consecutive_AT_scoring", "method", "mismatch_penalty_type", "fasta", "output_dir"):
            continue
        assert _params[key] >= 0.0, "Params must be a non-negative integer."

    params = Params(**_params)
    assert Path(fasta).expanduser().resolve().is_file(), f"No file {fasta} was found."
    
    print(colored("Process parameters", "magenta"))
    print(colored("*" * 25, "magenta"))
    for key, value in params.__new_dict__.items():
        print(colored(f"{key}: {value}", "magenta"))

    print(colored("*" * 25, "magenta"))
    logging.basicConfig(
        level=logging.WARNING,
        filemode="a",
        # filename="zdna_extractions.log",
        format="%(levelname)s:%(message)s"
    )
    transform(fasta, params=params)

if __name__ == "__main__":
    main()
