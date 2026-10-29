# ClinVar somatic/oncogenicity classifications for TSO500 variants


import argparse
import glob
import os
import re
import sys
import time
from typing import Dict, List, Optional

import pandas as pd
import requests

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
BATCH_SIZE = 100

GENE_COL = "gene"
CLINVAR_ID_COL = "clinvar.sig"

DENOVO_DIR = "<SMALL_VARIANT_DIR>"
RECURRENT_DIR = "<RECURRENT_SMALL_VARIANT_DIR>"


def sample_id_from_filename(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    m = re.match(r"^([A-Za-z0-9]+)_", stem)
    return m.group(1) if m else stem


def scan_folder(d: str) -> List[str]:
    if not d:
        return []
    if not os.path.isdir(d):
        print(f"Warning: folder does not exist -- skipped: {d}", file=sys.stderr)
        return []
    found = []
    for ext in ("*.xlsx", "*.xls"):
        found += glob.glob(os.path.join(d, "**", ext), recursive=True)
    found = [f for f in found if not os.path.basename(f).startswith("~$")]
    return found


def collect_input_files(explicit_files: List[str], cohort_dirs: Dict[str, str], extra_dirs: List[str]):
    entries = []

    for f in explicit_files:
        entries.append((f, "Unknown"))

    for d, cohort in cohort_dirs.items():
        found = scan_folder(d)
        print(f"  [{cohort}] {len(found)} file(s) found in {d}")
        for f in found:
            entries.append((f, cohort))

    for d in extra_dirs:
        found = scan_folder(d)
        print(f"  [Unknown] {len(found)} file(s) found in {d}")
        for f in found:
            entries.append((f, "Unknown"))

    seen = set()
    out = []
    for f, cohort in entries:
        if f not in seen:
            seen.add(f)
            out.append((f, cohort))
    return out


def load_variant_rows(files_with_cohort):
    frames = []
    sample_index_rows = []

    for f, cohort in files_with_cohort:
        sample = sample_id_from_filename(f)
        sample_index_rows.append({
            "sample": sample,
            "cohort": cohort,
            "source_file": os.path.basename(f),
        })

        df = pd.read_excel(f)
        missing = [c for c in (GENE_COL, CLINVAR_ID_COL) if c not in df.columns]
        if missing:
            print(f"Warning: {os.path.basename(f)} is missing column(s) {missing} -- skipped", file=sys.stderr)
            continue

        sub = df[[GENE_COL, CLINVAR_ID_COL]].copy()
        sub = sub.dropna(subset=[CLINVAR_ID_COL])
        if sub.empty:
            continue

        sub["sample"] = sample
        sub["cohort"] = cohort
        sub["source_file"] = os.path.basename(f)
        sub["variation_id"] = sub[CLINVAR_ID_COL].astype(float).astype("Int64").astype(str)
        frames.append(sub[["sample", "cohort", GENE_COL, "variation_id", "source_file"]])

    sample_index = pd.DataFrame(sample_index_rows).drop_duplicates(subset=["sample"])

    if not frames:
        variant_df = pd.DataFrame(columns=["sample", "cohort", GENE_COL, "variation_id", "source_file"])
    else:
        variant_df = pd.concat(frames, ignore_index=True)

    return variant_df, sample_index


def fetch_somatic_classifications(variation_ids: List[str], api_key: Optional[str]) -> Dict[str, dict]:
    results: Dict[str, dict] = {}
    unique_ids = sorted(set(variation_ids), key=int)

    for i in range(0, len(unique_ids), BATCH_SIZE):
        batch = unique_ids[i : i + BATCH_SIZE]
        params = {"db": "clinvar", "id": ",".join(batch), "retmode": "json"}
        if api_key:
            params["api_key"] = api_key

        for attempt in range(3):
            try:
                resp = requests.get(ESUMMARY_URL, params=params, timeout=30)
                resp.raise_for_status()
                data = resp.json()
                break
            except requests.exceptions.RequestException as e:
                wait = 2 ** attempt
                print(f"Warning: batch {i // BATCH_SIZE + 1} failed ({e}); "
                      f"retrying in {wait}s ({attempt + 1}/3) ...", file=sys.stderr)
                time.sleep(wait)
        else:
            print(f"Warning: batch {i // BATCH_SIZE + 1} failed after 3 attempts -- "
                  f"skipping {len(batch)} ID(s): {batch}", file=sys.stderr)
            continue

        result_block = data.get("result", {})
        for uid in result_block.get("uids", []):
            results[uid] = result_block.get(uid, {})

        time.sleep(0.35 if not api_key else 0.11)

    return results


def extract_somatic_calling(record: dict) -> Optional[dict]:
    for key in ("oncogenicity_classification", "clinical_impact_classification"):
        block = record.get(key)
        if isinstance(block, dict) and block.get("description"):
            return {
                "somatic_calling": block.get("description", ""),
                "classification_type": "oncogenicity" if key == "oncogenicity_classification" else "clinical_impact",
                "review_status": block.get("review_status", ""),
            }
    return None


def get_condition_names(record: dict) -> str:
    trait_names = []
    for key in ("oncogenicity_classification", "clinical_impact_classification"):
        block = record.get(key)
        if isinstance(block, dict):
            for t in block.get("trait_set", []) or []:
                name = t.get("trait_name")
                if name:
                    trait_names.append(name)
    return "; ".join(dict.fromkeys(trait_names))


def main():
    parser = argparse.ArgumentParser(description="Somatic/oncogenicity ClinVar calls for TSO500 filtered files.")
    parser.add_argument("files", nargs="*", help="Explicit .xlsx file paths (cohort tagged 'Unknown')")
    parser.add_argument("--denovo-dir", default=DENOVO_DIR, help="De novo cohort folder (recursively scanned)")
    parser.add_argument("--recurrent-dir", default=RECURRENT_DIR, help="Recurrent cohort folder (recursively scanned)")
    parser.add_argument("--dir", nargs="*", default=[], help="Extra folder(s) to scan, cohort tagged 'Unknown'")
    parser.add_argument("--api-key", help="NCBI API key (optional)")
    parser.add_argument("--out", default="somatic_clinvar_results.csv")
    parser.add_argument("--samples-out", default="samples_index.csv")
    args = parser.parse_args()

    cohort_dirs = {}
    if args.denovo_dir:
        cohort_dirs[args.denovo_dir] = "De novo"
    if args.recurrent_dir:
        cohort_dirs[args.recurrent_dir] = "Recurrence"

    print("Scanning folders...")
    files_with_cohort = collect_input_files(args.files, cohort_dirs, args.dir)
    if not files_with_cohort:
        parser.error("No input files found. Check --denovo-dir / --recurrent-dir, or pass file paths / --dir folders.")

    print(f"Reading {len(files_with_cohort)} file(s) total...")
    variant_df, sample_index = load_variant_rows(files_with_cohort)

    sample_index.to_csv(args.samples_out, index=False)
    print(f"Saved {len(sample_index)} sample(s) (all cohorts, incl. zero-hit patients) to {args.samples_out}")
    print(sample_index["cohort"].value_counts().rename("n_files").to_string())

    if variant_df.empty:
        print("No rows with a clinvar.sig ID were found across the given files.")
        pd.DataFrame(columns=[
            "sample", "cohort", "gene", "clinvar_variation_id", "somatic_calling",
            "classification_type", "review_status", "cancer_type_condition", "source_file",
        ]).to_csv(args.out, index=False)
        return

    print(f"Found {len(variant_df)} variant rows with a ClinVar ID "
          f"({variant_df['variation_id'].nunique()} unique variants). Querying ClinVar...")

    summaries = fetch_somatic_classifications(variant_df["variation_id"].tolist(), args.api_key)

    rows = []
    for _, r in variant_df.iterrows():
        record = summaries.get(r["variation_id"])
        if not record:
            continue
        somatic = extract_somatic_calling(record)
        if not somatic:
            continue
        rows.append({
            "sample": r["sample"],
            "cohort": r["cohort"],
            "gene": r[GENE_COL],
            "clinvar_variation_id": r["variation_id"],
            "somatic_calling": somatic["somatic_calling"],
            "classification_type": somatic["classification_type"],
            "review_status": somatic["review_status"],
            "cancer_type_condition": get_condition_names(record),
            "source_file": r["source_file"],
        })

    result_df = pd.DataFrame(rows)

    if result_df.empty:
        print("None of the queried variants had a somatic/oncogenicity classification in ClinVar.")
    else:
        print(result_df[["sample", "cohort", "gene", "somatic_calling"]].to_string(index=False))

    result_df.to_csv(args.out, index=False)
    print(f"\nSaved {len(result_df)} somatic-classified rows to {args.out}")


if __name__ == "__main__":
    main()
