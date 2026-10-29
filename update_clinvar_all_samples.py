# Re-query current ClinVar classifications for all TSO500 samples

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
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
BATCH_SIZE = 100

DENOVO_DIR = "<SMALL_VARIANT_DIR>"
RECURRENT_DIR = "<RECURRENT_SMALL_VARIANT_DIR>"

GENE_COL = "gene"
CHR_COL = "chr"
POS_COL = "POS"
REF_COL = "REF"
ALT_COL = "ALT"
AMINO_COL = "aminos.ann"
CDNA_COL = "cDNA.ann"
CLINVAR_ID_COL = "clinvar.sig"
OLD_SIG_COL = "clinvar.stat"
OLD_CONDITION_COL = "clinvar.CLNDN"


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


def collect_input_files(explicit_files: List[str], cohort_dirs: Dict[str, str]):
    entries = [(f, "Unknown") for f in explicit_files]

    for d, cohort in cohort_dirs.items():
        found = scan_folder(d)
        print(f"  [{cohort}] {len(found)} file(s) found in {d}")
        for f in found:
            entries.append((f, cohort))

    seen = set()
    out = []
    for f, cohort in entries:
        if f not in seen:
            seen.add(f)
            out.append((f, cohort))
    return out


def load_variant_rows(files_with_cohort) -> (pd.DataFrame, pd.DataFrame):
    frames = []
    sample_index_rows = []

    for f, cohort in files_with_cohort:
        sample = sample_id_from_filename(f)
        sample_index_rows.append({"sample": sample, "cohort": cohort, "source_file": os.path.basename(f)})

        try:
            df = pd.read_excel(f)
        except Exception as e:
            print(f"Warning: could not read {os.path.basename(f)} ({e}) -- skipped", file=sys.stderr)
            continue

        missing = [c for c in (GENE_COL,) if c not in df.columns]
        if missing:
            print(f"Warning: {os.path.basename(f)} is missing column(s) {missing} -- skipped", file=sys.stderr)
            continue

        keep_cols = [c for c in (
            GENE_COL, CHR_COL, POS_COL, REF_COL, ALT_COL,
            AMINO_COL, CDNA_COL, CLINVAR_ID_COL, OLD_SIG_COL, OLD_CONDITION_COL,
        ) if c in df.columns]

        sub = df[keep_cols].copy()
        if sub.empty:
            continue

        sub["sample"] = sample
        sub["cohort"] = cohort
        sub["source_file"] = os.path.basename(f)

        if CLINVAR_ID_COL in sub.columns:
            has_id = sub[CLINVAR_ID_COL].notna()
        else:
            sub[CLINVAR_ID_COL] = pd.NA
            has_id = pd.Series(False, index=sub.index)

        sub["variation_id"] = pd.NA
        sub.loc[has_id, "variation_id"] = (
            sub.loc[has_id, CLINVAR_ID_COL].astype(float).astype("Int64").astype(str)
        )
        sub["id_source"] = ""
        sub.loc[has_id, "id_source"] = "clinvar.sig column"

        frames.append(sub)

    sample_index = pd.DataFrame(sample_index_rows).drop_duplicates(subset=["sample"])
    variant_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return variant_df, sample_index


def find_variation_id_by_hgvs(hgvs: str, api_key: Optional[str]) -> dict:
    params = {"db": "clinvar", "term": hgvs, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key

    for attempt in range(3):
        try:
            resp = requests.get(ESEARCH_URL, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            break
        except requests.exceptions.RequestException as e:
            wait = 2 ** attempt
            print(f"Warning: esearch for '{hgvs}' failed ({e}); retrying in {wait}s ...", file=sys.stderr)
            time.sleep(wait)
    else:
        return {"variation_id": None, "n_matches": -1, "note": "esearch failed after 3 attempts"}

    idlist = data.get("esearchresult", {}).get("idlist", [])
    if len(idlist) == 0:
        return {"variation_id": None, "n_matches": 0, "note": "no ClinVar match for this HGVS string"}
    if len(idlist) == 1:
        return {"variation_id": idlist[0], "n_matches": 1, "note": ""}
    return {"variation_id": None, "n_matches": len(idlist),
            "note": f"ambiguous -- {len(idlist)} candidate IDs: {', '.join(idlist)}"}


def fill_missing_ids_via_hgvs(variants: pd.DataFrame, api_key: Optional[str]) -> pd.DataFrame:
    variants["hgvs_lookup_note"] = ""
    if variants.empty:
        return variants

    missing_mask = variants["variation_id"].isna()
    if CDNA_COL not in variants.columns:
        if missing_mask.any():
            print("No 'cDNA.ann' column available -- cannot attempt HGVS fallback.")
        return variants

    n_to_try = int(missing_mask.sum())
    if n_to_try == 0:
        return variants

    print(f"Attempting HGVS-based esearch lookup for {n_to_try} row(s) with no VariationID "
          f"(this is the slow step -- one request per row, rate-limited) ...")

    done = 0
    for idx in variants[missing_mask].index:
        hgvs = variants.at[idx, CDNA_COL] if CDNA_COL in variants.columns else None
        if pd.isna(hgvs) or not str(hgvs).strip():
            variants.at[idx, "hgvs_lookup_note"] = "no HGVS string available to search"
            continue

        hgvs = str(hgvs).strip()
        result = find_variation_id_by_hgvs(hgvs, api_key)

        if result["variation_id"]:
            variants.at[idx, "variation_id"] = result["variation_id"]
            variants.at[idx, "id_source"] = f"HGVS esearch ({hgvs})"
        else:
            variants.at[idx, "hgvs_lookup_note"] = result["note"]

        done += 1
        if done % 25 == 0:
            print(f"  ... {done}/{n_to_try} HGVS lookups done")

        time.sleep(0.35 if not api_key else 0.11)

    return variants


def fetch_current_classifications(variation_ids, api_key: Optional[str]) -> Dict[str, dict]:
    results: Dict[str, dict] = {}
    unique_ids = sorted(set(variation_ids), key=int)
    print(f"Querying esummary for {len(unique_ids)} unique VariationID(s) across all samples "
          f"({(len(unique_ids) + BATCH_SIZE - 1) // BATCH_SIZE} batch(es)) ...")

    for i in range(0, len(unique_ids), BATCH_SIZE):
        batch = unique_ids[i:i + BATCH_SIZE]
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
                print(f"Warning: batch {i // BATCH_SIZE + 1} failed ({e}); retrying in {wait}s ({attempt + 1}/3) ...",
                      file=sys.stderr)
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


def get_germline(record: dict) -> dict:
    block = record.get("germline_classification")
    if isinstance(block, dict):
        return {
            "current_germline_sig": block.get("description", ""),
            "current_germline_review_status": block.get("review_status", ""),
        }
    return {"current_germline_sig": "", "current_germline_review_status": ""}


def get_somatic(record: dict) -> dict:
    for key, label in (
        ("oncogenicity_classification", "oncogenicity"),
        ("clinical_impact_classification", "clinical_impact"),
    ):
        block = record.get(key)
        if isinstance(block, dict) and block.get("description"):
            traits = []
            for t in block.get("trait_set", []) or []:
                name = t.get("trait_name")
                if name:
                    traits.append(name)
            return {
                "current_somatic_calling": block.get("description", ""),
                "somatic_classification_type": label,
                "current_somatic_review_status": block.get("review_status", ""),
                "current_somatic_condition": "; ".join(dict.fromkeys(traits)),
            }
    return {
        "current_somatic_calling": "",
        "somatic_classification_type": "",
        "current_somatic_review_status": "",
        "current_somatic_condition": "",
    }


def main():
    parser = argparse.ArgumentParser(description="Re-query current ClinVar classifications for all TSO500 filtered files in both cohorts.")
    parser.add_argument("files", nargs="*", help="Explicit .xlsx file paths (cohort tagged 'Unknown'), in addition to the two folders")
    parser.add_argument("--denovo-dir", default=DENOVO_DIR, help="De novo cohort folder (recursively scanned)")
    parser.add_argument("--recurrent-dir", default=RECURRENT_DIR, help="Recurrent cohort folder (recursively scanned)")
    parser.add_argument("--api-key", help="NCBI API key (optional, raises rate limit 3/s -> 10/s)")
    parser.add_argument("--no-hgvs-fallback", action="store_true", help="Skip the esearch-by-HGVS fallback step")
    parser.add_argument("--out", default="clinvar_updated_all_samples.csv")
    parser.add_argument("--samples-out", default="samples_index.csv")
    args = parser.parse_args()

    cohort_dirs = {}
    if args.denovo_dir:
        cohort_dirs[args.denovo_dir] = "De novo"
    if args.recurrent_dir:
        cohort_dirs[args.recurrent_dir] = "Recurrence"

    print("Scanning folders...")
    files_with_cohort = collect_input_files(args.files, cohort_dirs)
    if not files_with_cohort:
        parser.error("No input files found. Check --denovo-dir / --recurrent-dir, or pass file paths directly.")

    print(f"Reading {len(files_with_cohort)} file(s) total...")
    variants, sample_index = load_variant_rows(files_with_cohort)

    sample_index.to_csv(args.samples_out, index=False)
    print(f"Saved {len(sample_index)} sample(s) (all cohorts) to {args.samples_out}")
    if not sample_index.empty:
        print(sample_index["cohort"].value_counts().rename("n_files").to_string())

    if variants.empty:
        print("No variant rows were found across the given files.")
        return

    n_total = len(variants)
    n_no_id = int(variants["variation_id"].isna().sum())
    print(f"\n{n_total} variant row(s) total across all samples; "
          f"{n_total - n_no_id} already have a VariationID, {n_no_id} do not.")

    if not args.no_hgvs_fallback:
        variants = fill_missing_ids_via_hgvs(variants, args.api_key)
    else:
        variants["hgvs_lookup_note"] = ""

    queryable = variants.dropna(subset=["variation_id"]).copy()
    still_missing = len(variants) - len(queryable)
    print(f"\n{len(queryable)} of {len(variants)} row(s) have a ClinVar VariationID "
          f"(directly or via HGVS lookup) and will be queried against live ClinVar.")
    if still_missing:
        print(f"{still_missing} row(s) remain without an ID -- see 'hgvs_lookup_note' in the output CSV.")

    summaries = {}
    if not queryable.empty:
        summaries = fetch_current_classifications(queryable["variation_id"].tolist(), args.api_key)

    rows = []
    for _, r in queryable.iterrows():
        record = summaries.get(r["variation_id"])
        row = {
            "sample": r.get("sample", ""),
            "cohort": r.get("cohort", ""),
            "gene": r.get(GENE_COL, ""),
            "chr": r.get(CHR_COL, ""),
            "POS": r.get(POS_COL, ""),
            "REF": r.get(REF_COL, ""),
            "ALT": r.get(ALT_COL, ""),
            "aminos.ann": r.get(AMINO_COL, ""),
            "cDNA.ann": r.get(CDNA_COL, ""),
            "clinvar_variation_id": r["variation_id"],
            "id_source": r.get("id_source", ""),
            "old_clinvar_stat": r.get(OLD_SIG_COL, ""),
            "old_clinvar_condition": r.get(OLD_CONDITION_COL, ""),
            "source_file": r.get("source_file", ""),
        }
        if not record:
            row.update({
                "current_germline_sig": "NOT FOUND",
                "current_germline_review_status": "",
                "current_somatic_calling": "",
                "somatic_classification_type": "",
                "current_somatic_review_status": "",
                "current_somatic_condition": "",
                "hgvs_lookup_note": "",
            })
        else:
            row.update(get_germline(record))
            row.update(get_somatic(record))
            row["hgvs_lookup_note"] = ""
        rows.append(row)

    result_df = pd.DataFrame(rows)

    unmatched = variants[variants["variation_id"].isna()].copy()
    if not unmatched.empty:
        unmatched_rows = []
        for _, r in unmatched.iterrows():
            unmatched_rows.append({
                "sample": r.get("sample", ""),
                "cohort": r.get("cohort", ""),
                "gene": r.get(GENE_COL, ""),
                "chr": r.get(CHR_COL, ""),
                "POS": r.get(POS_COL, ""),
                "REF": r.get(REF_COL, ""),
                "ALT": r.get(ALT_COL, ""),
                "aminos.ann": r.get(AMINO_COL, ""),
                "cDNA.ann": r.get(CDNA_COL, ""),
                "clinvar_variation_id": "",
                "id_source": "NOT FOUND (checked clinvar.sig + HGVS esearch)",
                "old_clinvar_stat": r.get(OLD_SIG_COL, ""),
                "old_clinvar_condition": r.get(OLD_CONDITION_COL, ""),
                "source_file": r.get("source_file", ""),
                "current_germline_sig": "",
                "current_germline_review_status": "",
                "current_somatic_calling": "",
                "somatic_classification_type": "",
                "current_somatic_review_status": "",
                "current_somatic_condition": "",
                "hgvs_lookup_note": r.get("hgvs_lookup_note", ""),
            })
        result_df = pd.concat([result_df, pd.DataFrame(unmatched_rows)], ignore_index=True)

    result_df["germline_changed"] = (
        result_df["old_clinvar_stat"].fillna("").astype(str).str.replace(" ", "").str.lower()
        != result_df["current_germline_sig"].fillna("").astype(str).str.replace(" ", "").str.lower()
    )

    result_df.to_csv(args.out, index=False)
    print(f"\nSaved {len(result_df)} row(s) across {result_df['sample'].nunique()} sample(s) to {args.out}")

    n_changed = int(result_df["germline_changed"].sum())
    if n_changed:
        print(f"{n_changed} variant(s) show a DIFFERENT germline classification today vs. the original file -- review these.")
    n_somatic = int((result_df["current_somatic_calling"] != "").sum())
    print(f"{n_somatic} variant(s) now have a somatic/oncogenicity classification in ClinVar "
          f"(not present in the original germline-only extraction).")

    if not result_df.empty:
        print("\nBy cohort:")
        print(result_df.groupby("cohort")["sample"].nunique().rename("n_samples_with_matched_variants").to_string())


if __name__ == "__main__":
    main()
