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
import bisect

@define(kw_only=True, slots=True, frozen=True)
class ZDNASubmissionForm:
    recordSeq: str = field(converter=str)
    recordID: str = field(converter=str)
    input_fasta: os.PathLike[str] = field(converter=str)

def timeit(func: Callable) -> Callable:
    """
    Simple decorator to measure execution time for a function.
    """
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
    """
    Yields (record_id, record_sequence).
    Handles gzipped FASTA if extension ends with .gz
    """
    if Path(fasta).name.endswith(".gz"):
        file = gzip.open(fasta, 'rt')
    else:
        file = open(fasta, encoding='utf-8', mode='r')

    for record in SimpleFastaParser(file):
        # record[0] -> header line
        # record[1] -> sequence
        yield record[0].split(' ')[0], record[1]
    file.close()

def _extract(ID: int, sequence: str, params: Params) -> pd.DataFrame:
    """
    For a given sequence, calculates Z-DNA regions above threshold and returns
    them as a DataFrame.
    """
    zdna_string = ZDNACalculatorSeq(data=sequence, params=params)
    subarrays_detected = zdna_string.subarrays_above_threshold()
    subarrays_detected = pd.DataFrame(subarrays_detected, columns=params.headers[1:-1])
    scoring_array = zdna_string.scoring_array

    if not len(subarrays_detected):
        # If no subarrays found, produce a row of NaNs for consistency
        null_subarrays = pd.DataFrame([[np.nan]*len(params.headers)], columns=params.headers)
        null_subarrays['Chromosome'] = null_subarrays['Chromosome'].astype(object)
        null_subarrays.loc[:, "Chromosome"] = ID
        null_subarrays.loc[:, "totalSequenceScore"] = sum(scoring_array)
        return null_subarrays

    subarrays_detected.loc[:, "Chromosome"] = ID
    subarrays_detected.loc[:, "totalSequenceScore"] = sum(scoring_array)
    return subarrays_detected[params.headers]

def _subextraction(submission_forms: list[ZDNASubmissionForm], params: Params) -> pd.DataFrame:
    """
    Processes a chunk of sequences (submission_forms) in a single process,
    returning a DataFrame with the combined Z-DNA results for those sequences.
    """
    outputs = []
    for form in submission_forms:
        try:
            present_df = _extract(ID=form.recordID, sequence=form.recordSeq, params=params)
            outputs.append(present_df)
        except IndexError as e:
            logging.error(
                f"Failed to process input fasta '{Path(form.input_fasta).name}' "
                f"on chromosome {str(form.recordID)}. "
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
    """
    Splits 'tasks' into 'total_buckets' roughly even parts
    for parallel processing.
    """
    then = time.perf_counter()
    total = len(tasks)
    if total_buckets <= 0:
        total_buckets = 1
    if total_buckets > total:
        # If we have more buckets than tasks, some buckets will be empty
        total_buckets = total

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
    """
    Main entry point for extracting Z-DNA subarrays from a FASTA file
    using parallel processing according to params.n_jobs.
    """
    then = time.perf_counter()
    zdna_df = pd.DataFrame(columns=params.headers)
    outputs: list[pd.DataFrame] = []
    n_jobs = params.n_jobs
    assert isinstance(n_jobs, int) and n_jobs > 0, "Number of jobs must be a positive int."

    submission_forms: list[ZDNASubmissionForm] = []

    # Gather sequences
    for rec_id, rec_seq in read_fasta(fasta):
        submission_forms.append(
            ZDNASubmissionForm(
                recordSeq=rec_seq.upper(),
                recordID=rec_id,
                input_fasta=fasta
            )
        )

    total_submission_forms = len(submission_forms)

    # Assign tasks for parallel extraction
    assigned_tasks = assign_tasks(submission_forms, n_jobs)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=n_jobs, mp_context=mp.get_context("spawn")
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

    # If display_sequence_score == 0, drop totalSequenceScore and also drop rows of all NaN
    if params.display_sequence_score == 0:
        zdna_df.drop(columns=['totalSequenceScore'], inplace=True, errors='ignore')
        zdna_df.dropna(axis=0, how='all', inplace=True)

    return zdna_df

def parse_consecutive_AT_scoring(value: str) -> tuple[float, ...]:
    """
    Parse a comma-separated string of floats (like "3.0, 1.5, 0.7") into a tuple[float,...].
    """
    try:
        values = [float(x.strip()) for x in value.split(',')]
        return tuple(values)
    except ValueError as e:
        raise argparse.ArgumentTypeError(f"Invalid consecutive_AT_scoring format: {e}")

def parse_gff_file(gff_path: Path) -> pd.DataFrame:
    """
    Reads a GFF3 file, filters for lines where the 3rd column == 'gene'.
    Returns a DataFrame with columns:
      Chromosome, gene_start, gene_end, strand, gene_id, gene_biotype, locus_tag

    Important: GFF is usually 1-based inclusive. We convert to 0-based inclusive
    by subtracting 1 from both start and end so that it matches
    the 0-based indexing used in the Z-DNA pipeline.
    """
    rows = []
    with open(gff_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            if feature_type.lower() != 'gene':
                # skip everything but actual 'gene' features
                continue

            chrom = parts[0]
            # Convert 1-based [start, end] to 0-based inclusive
            start = int(parts[3]) - 1
            end = int(parts[4]) - 1
            strand = parts[6]

            # The last column has attributes in the form: key1=value1;key2=value2;...
            attributes = parts[8].split(';')
            attr_dict = {}
            for kv in attributes:
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attr_dict[k.strip()] = v.strip()

            gene_id = attr_dict.get('ID', 'NA')
            gene_biotype = attr_dict.get('gene_biotype', 'NA')
            locus_tag = attr_dict.get('locus_tag', gene_id)

            rows.append({
                "Chromosome": chrom,
                "gene_start": start,
                "gene_end": end,
                "strand": strand,
                "gene_id": gene_id,
                "gene_biotype": gene_biotype,
                "locus_tag": locus_tag
            })

    df = pd.DataFrame(rows)
    # Sort by gene_start ascending (this helps us do binary search quickly)
    df.sort_values(by=["Chromosome", "gene_start"], inplace=True)
    return df.reset_index(drop=True)

def compute_distances(ir_start: int,
                      ir_end: int,
                      gene_start: int,
                      gene_end: int,
                      strand: str) -> tuple[int, int, int]:
    """
    Returns (distance, distance_from_TSS, distance_from_TES) for an IR (ir_start, ir_end)
    and a gene (gene_start, gene_end, strand).

    - For a plus-strand gene, TSS = gene_start, TES = gene_end.
    - For a minus-strand gene, TSS = gene_end,   TES = gene_start.
    - distance=0 if there's any overlap; otherwise it's the minimal gap
      between the IR and the gene.
    """
    # Check overlap: IR_end < gene_start means IR is entirely to the left,
    # or IR_start > gene_end means IR is entirely to the right.
    overlap = not (ir_end < gene_start or ir_start > gene_end)

    if strand == '+':
        TSS = gene_start
        TES = gene_end
    else:  # '-'
        TSS = gene_end
        TES = gene_start

    # distance_from_TSS
    if strand == '+':
        if overlap:
            dist_tss = max(ir_start - TSS, 0)
        else:
            # IR fully to the left
            if ir_end < gene_start:
                dist_tss = gene_start - ir_end
            else:
                # IR fully to the right
                dist_tss = ir_start - TSS
    else:  # '-'
        if overlap:
            dist_tss = max(TSS - ir_end, 0) if TSS > ir_end else 0
        else:
            if ir_end < gene_start:
                dist_tss = TSS - ir_end if TSS > ir_end else 0
            else:
                dist_tss = ir_start - TSS if ir_start > TSS else 0

    # distance_from_TES
    if strand == '+':
        if overlap:
            dist_tes = max(gene_end - ir_end, 0)
        else:
            if ir_start > gene_end:
                dist_tes = ir_start - gene_end
            else:
                dist_tes = gene_end - ir_end
    else:  # '-'
        if overlap:
            dist_tes = max(ir_start - TES, 0) if ir_start > TES else 0
        else:
            if ir_end < gene_start:
                dist_tes = gene_start - ir_end
            else:
                dist_tes = ir_start - gene_start

    # Overall distance (0 if overlap)
    if overlap:
        distance = 0
    else:
        # Minimal gap from gene boundaries
        left_gap = gene_start - ir_end if ir_end < gene_start else None
        right_gap = ir_start - gene_end if ir_start > gene_end else None
        possible_gaps = [abs(g) for g in (left_gap, right_gap) if g is not None]
        distance = min(possible_gaps) if possible_gaps else 0

    return (distance, dist_tss, dist_tes)

def annotate_closest_gene_in_chromosome(chrom_ir_df: pd.DataFrame,
                                        chrom_gene_df: pd.DataFrame) -> pd.DataFrame:
    """
    For IRs in a single chromosome, find the closest gene from chrom_gene_df (also for that chromosome).
    Both DataFrames are assumed sorted by 'Start' (IR) and 'gene_start' (genes) respectively.
    Uses binary search to find candidate genes.
    Returns a new DataFrame with annotation columns appended:
        gene_start, gene_end, gene_id, gene_biotype, strand,
        distance, distance_from_TSS, distance_from_TES
    """
    if chrom_ir_df.empty or chrom_gene_df.empty:
        # Just fill with NaNs
        chrom_ir_df["gene_start"] = np.nan
        chrom_ir_df["gene_end"] = np.nan
        chrom_ir_df["gene_id"] = np.nan
        chrom_ir_df["gene_biotype"] = np.nan
        chrom_ir_df["strand"] = np.nan
        chrom_ir_df["distance"] = np.nan
        chrom_ir_df["distance_from_TSS"] = np.nan
        chrom_ir_df["distance_from_TES"] = np.nan
        return chrom_ir_df

    # Sort IR data by IR start
    chrom_ir_df = chrom_ir_df.sort_values(by="Start").reset_index(drop=True)

    # We'll store gene_start in a list for binary search
    gene_starts = chrom_gene_df["gene_start"].tolist()

    new_cols = {
        "gene_start": [],
        "gene_end": [],
        "gene_id": [],
        "gene_biotype": [],
        "strand": [],
        "distance": [],
        "distance_from_TSS": [],
        "distance_from_TES": []
    }

    for idx, row in chrom_ir_df.iterrows():
        ir_start = int(row["Start"])
        ir_end = int(row["End"])

        # Binary-search the position of ir_start in gene_starts
        i = bisect.bisect_left(gene_starts, ir_start)

        best_dist = None
        best_entry = None
        best_dist_tss = None
        best_dist_tes = None

        # We'll check the gene at index i, plus neighbors i-1, i+1,
        # because the closest gene might be just before or after the insertion point.
        candidates_idx = []
        for offset in [i, i-1, i+1]:
            if 0 <= offset < len(chrom_gene_df):
                candidates_idx.append(offset)

        for cidx in candidates_idx:
            g = chrom_gene_df.iloc[cidx]
            distance, dist_tss, dist_tes = compute_distances(
                ir_start, ir_end,
                g["gene_start"], g["gene_end"],
                g["strand"]
            )
            if (best_dist is None) or (distance < best_dist):
                best_dist = distance
                best_dist_tss = dist_tss
                best_dist_tes = dist_tes
                best_entry = g

        if best_entry is not None:
            new_cols["gene_start"].append(best_entry["gene_start"])
            new_cols["gene_end"].append(best_entry["gene_end"])
            new_cols["gene_id"].append(best_entry["gene_id"])
            new_cols["gene_biotype"].append(best_entry["gene_biotype"])
            new_cols["strand"].append(best_entry["strand"])
            new_cols["distance"].append(best_dist)
            new_cols["distance_from_TSS"].append(best_dist_tss)
            new_cols["distance_from_TES"].append(best_dist_tes)
        else:
            new_cols["gene_start"].append(np.nan)
            new_cols["gene_end"].append(np.nan)
            new_cols["gene_id"].append(np.nan)
            new_cols["gene_biotype"].append(np.nan)
            new_cols["strand"].append(np.nan)
            new_cols["distance"].append(np.nan)
            new_cols["distance_from_TSS"].append(np.nan)
            new_cols["distance_from_TES"].append(np.nan)

    for col, values in new_cols.items():
        chrom_ir_df[col] = values

    return chrom_ir_df

def process_chromosome_group(chrom_group, gene_groups):
    """
    Process a single (chrom_name, chrom_ir_df) group by finding
    the closest genes in gene_groups for that same chrom_name.
    Returns a (possibly empty) DataFrame with new annotation columns.
    """
    chrom_name, chrom_ir_df = chrom_group
    if chrom_name not in gene_groups:
        # No genes on this chromosome
        chrom_ir_df["gene_start"] = np.nan
        chrom_ir_df["gene_end"] = np.nan
        chrom_ir_df["gene_id"] = np.nan
        chrom_ir_df["gene_biotype"] = np.nan
        chrom_ir_df["strand"] = np.nan
        chrom_ir_df["distance"] = np.nan
        chrom_ir_df["distance_from_TSS"] = np.nan
        chrom_ir_df["distance_from_TES"] = np.nan
        return chrom_ir_df

    chrom_gene_df = gene_groups[chrom_name]
    return annotate_closest_gene_in_chromosome(chrom_ir_df, chrom_gene_df)

def process_chunk_of_chrom_groups(chunk, gene_groups):
    """
    For a list of chromosome groups (chunk), process each group
    and return a single DataFrame with all results concatenated.

    If the chunk is empty, returns an empty DataFrame.
    If all DataFrames are empty, also returns an empty DataFrame.
    """
    if not chunk:
        return pd.DataFrame()

    dfs = []
    for g in chunk:
        df_annotated = process_chromosome_group(g, gene_groups)
        if df_annotated is not None:
            dfs.append(df_annotated)

    # Filter out any empty dataframes before concatenating
    dfs_non_empty = [d for d in dfs if not d.empty]

    if not dfs_non_empty:
        return pd.DataFrame()
    return pd.concat(dfs_non_empty, ignore_index=True)

def annotate_closest_gene(zdna_df: pd.DataFrame, gff_df: pd.DataFrame, n_jobs: int) -> pd.DataFrame:
    """
    For each chromosome in zdna_df, find the closest gene from gff_df
    (filtered by that chromosome). Uses binary search on sorted gene
    coordinates for efficiency. Parallelizes at the chromosome level
    using n_jobs.
    """
    if zdna_df.empty or gff_df.empty:
        # Edge case: no data or no genes
        for col in [
            "gene_start","gene_end","gene_id","gene_biotype","strand",
            "distance","distance_from_TSS","distance_from_TES"
        ]:
            zdna_df[col] = np.nan
        return zdna_df

    # Group IR by chromosome
    zdna_groups = list(zdna_df.groupby("Chromosome"))
    # Group genes by chromosome
    gene_groups = dict(tuple(gff_df.groupby("Chromosome")))

    # Split the chromosome groups into chunks for parallel processing
    assigned_tasks = assign_tasks(zdna_groups, n_jobs)

    results = []
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=n_jobs, mp_context=mp.get_context("spawn")
    ) as executor:
        future_to_chunk = {
            executor.submit(process_chunk_of_chrom_groups, chunk, gene_groups): chunk
            for chunk in assigned_tasks
        }

        for future in concurrent.futures.as_completed(future_to_chunk):
            try:
                res = future.result()
                if res is not None and not res.empty:
                    results.append(res)
            except Exception as exc:
                logging.error(f"Chunk annotation failed with exception: {exc}")
                # Optionally re-raise if you want to halt execution
                # raise

    if results:
        annotated_df = pd.concat(results, ignore_index=True)
    else:
        # If all chunks were empty or had no results
        annotated_df = pd.DataFrame(
            columns=zdna_df.columns.tolist() + [
                "gene_start","gene_end","gene_id","gene_biotype","strand",
                "distance","distance_from_TSS","distance_from_TES"
            ]
        )
    return annotated_df

@timeit
def transform(path: Path | str, params: Params, gff_file: Path | None = None) -> pd.DataFrame:
    """
    Loads a FASTA, calculates Z-DNA features, optionally annotates with gene info (if gff_file is provided),
    then saves the final CSV.
    """
    logging.info(f"Processing file '{Path(path).name}'...")
    output_dir = params.output_dir
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    # 1) Extract Z-DNA
    zdna_df = extract_zdna_v2(path, params)

    # 2) If a GFF file is specified, parse and annotate
    if gff_file is not None and gff_file.is_file() and not zdna_df.empty:
        print(colored(f"Annotating Z-DNA results with genes from '{gff_file.name}'...", "cyan"))
        gff_df = parse_gff_file(gff_file)
        zdna_df = annotate_closest_gene(zdna_df, gff_df, params.n_jobs)

    # 3) Save final CSV
    print(colored("Saving Z-DNA dataframe...", "magenta"))
    output_file = output_dir.joinpath(f"{Path(path).stem}_zdna_score.csv")
    zdna_df.to_csv(output_file, mode="w", index=False)
    logging.info(f"File '{Path(path).name}' has been processed successfully.")
    return zdna_df

def main():
    parser = argparse.ArgumentParser(
        description="""Given a FASTA file and parameters, calculates ZDNA for each sequence present.
                       Optionally, supply --gff_file to annotate hits with nearest genes."""
    )
    parser.add_argument("--fasta", type=str, default="test_file.fna", help="Path to the FASTA file.")
    parser.add_argument("--GC_weight", type=float, default=Params.GC_weight,
                        help=f"Weight given to GC and CG transitions. Default = {Params.GC_weight}")
    parser.add_argument("--AT_weight", type=float, default=Params.AT_weight,
                        help=f"Weight given to AT and TA transitions. Default = {Params.AT_weight}")
    parser.add_argument("--GT_weight", type=float, default=Params.GT_weight,
                        help=f"Weight given to GT and TG transitions. Default = {Params.GT_weight}")
    parser.add_argument("--AC_weight", type=float, default=Params.AC_weight,
                        help=f"Weight given to AC and CA transitions. Default = {Params.AC_weight}")
    parser.add_argument("--mismatch_penalty_starting_value", type=int,
                        default=Params.mismatch_penalty_starting_value,
                        help="Penalty for the first non purine/pyrimidine transition encountered.")
    parser.add_argument("--mismatch_penalty_linear_delta", type=int,
                        default=Params.mismatch_penalty_linear_delta,
                        help="Rate of penalty increase for each subsequent mismatch if penalty type is linear.")
    parser.add_argument("--mismatch_penalty_type", choices=Params.mismatch_penalty_choices,
                        default=Params.mismatch_penalty_type,
                        help="Method of scaling the penalty for contiguous non purine/pyrimidine transitions.")
    parser.add_argument("--n_jobs", type=int, default=Params.n_jobs,
                        help="Number of processes to use. Default = -1 (all available cores).")
    parser.add_argument("--threshold", type=int, default=Params.threshold,
                        help="Scoring threshold above which a subarray is considered Z-DNA forming.")
    parser.add_argument("--consecutive_AT_scoring", type=parse_consecutive_AT_scoring,
                        default=Params.consecutive_AT_scoring,
                        help="Penalty array for consecutive AT repeats. Format like '3.0,2.0,1.0'.")
    parser.add_argument("--display_sequence_score", type=int, choices=[0, 1], default=0,
                        help="If 1, keep the totalSequenceScore column in the output.")
    parser.add_argument("--output_dir", type=str, default="zdna_extractions",
                        help="Output directory for CSV files.")
    parser.add_argument("--gff_file", type=str, default=None,
                        help="Optional GFF file for gene annotation. Only 'gene' features are used.")

    args = parser.parse_args()
    _params = vars(args)
    fasta = _params.pop("fasta")

    # Grab the gff_file path (may be None)
    gff_path = _params.pop("gff_file", None)
    gff_path = Path(gff_path) if gff_path else None

    # Basic param validation
    for key in _params:
        if key in ("consecutive_AT_scoring", "method", "mismatch_penalty_type", "output_dir"):
            continue
        assert _params[key] is not None, f"Parameter '{key}' cannot be None."
        if isinstance(_params[key], (int, float)) and key not in ("mismatch_penalty_type",):
            assert _params[key] >= 0.0, f"Param {key} must be a non-negative value."

    params = Params(**_params)
    fasta_path = Path(fasta).expanduser().resolve()
    assert fasta_path.is_file(), f"No file {fasta} was found."

    print(colored("Process parameters", "magenta"))
    print(colored("*" * 25, "magenta"))
    for key, value in params.__new_dict__.items():
        print(colored(f"{key}: {value}", "magenta"))
    print(colored("*" * 25, "magenta"))

    logging.basicConfig(
        level=logging.WARNING,
        filemode="a",
        format="%(levelname)s:%(message)s"
    )

    # Run main transform, optionally pass gff_path
    transform(path=fasta_path, params=params, gff_file=gff_path)

if __name__ == "__main__":
    main()
