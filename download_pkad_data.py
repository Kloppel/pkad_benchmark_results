#!/usr/bin/env python3
"""download_pkad_data.py — Read the three PKAD tables from the repo root, normalise
them to a common CSV schema, download crystal PDB structures via Biopython, and
download AlphaFold DB predictions for structures with a UniProt ID (PKAD3).

Source files (repo root):
  PKAD2_DOWNLOAD.xlsx_0.ods  — sheet "Wild Type", 1 461 rows
  PKAD3.csv                  — 1 805 rows (type column: single-letter D/E/H/K/C/Y/R)
  PKAD-R.csv                 — 1 025 rows (ResName already 3-letter)

Output layout (relative to this script):
  pkad_data/
    clean/
      PKAD2_clean.csv
      PKAD3_clean.csv
      PKAD-R_clean.csv
    super_table.csv
    pdbs/heavy_atom/<pdb_id>.pdb
    af3_pdbs/AF-<uniprot_id>-F1-model_v<N>.pdb
    missing_af3.csv

Usage
-----
    python3 download_pkad_data.py [--pdbs-only] [--no-af3] [--skip-existing]
                                  [--data-dir DIR]

Options
-------
  --pdbs-only      Skip table processing; only (re-)download PDB structures.
  --no-af3         Skip AlphaFold DB downloads.
  --skip-existing  Do not re-download files that already exist (default: on).
  --no-skip-existing  Force re-download everything.
  --data-dir DIR   Override pkad_data/ output directory.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HERE = Path(__file__).parent
REPO_ROOT = HERE.parents[1]          # benchmarking/pkad_benchmark → benchmarking → Karlsberg3/

PKAD2_SRC = REPO_ROOT / "PKAD2_DOWNLOAD.xlsx_0.ods"
PKAD3_SRC = REPO_ROOT / "PKAD3.csv"
PKADР_SRC = REPO_ROOT / "PKAD-R.csv"

# Common output schema column order
_SCHEMA = [
    "source_dataset", "pdb_id", "chain", "resname", "resid",
    "expt_pka", "expt_uncertainty",
    "uniprot_id", "uniprot_resid",
    "expt_method", "reference", "notes",
]

# Single-letter type → 3-letter residue name (PKAD3)
AA1_TO_3: dict[str, str] = {
    "D": "ASP", "E": "GLU", "H": "HIS", "K": "LYS",
    "C": "CYS", "Y": "TYR", "R": "ARG",
    "N": "NTE", "T": "CTE",
}

# ---------------------------------------------------------------------------
# Biopython imports (lazy — only needed when downloading crystal PDBs)
# ---------------------------------------------------------------------------

def _import_biopython():
    try:
        from Bio.PDB import PDBList, PDBParser, PDBIO, Select
        return PDBList, PDBParser, PDBIO, Select
    except ImportError:
        print("[ERROR] Biopython not installed. Run: pip install biopython", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# 1. Table readers
# ---------------------------------------------------------------------------

def read_pkad2_ods(path: Path) -> pd.DataFrame:
    """Read PKAD2 ODS (sheet 'Wild Type') and return normalised DataFrame."""
    print(f"[PKAD2] Reading {path.name} …")
    df = pd.read_excel(path, sheet_name="Wild Type", engine="odf", dtype=str)
    df = df.rename(columns={
        "PDB ID":                "pdb_id",
        "Res. Name":             "resname",
        "Chain":                 "chain",
        "Res. ID":               "resid",
        "Expt. pKa":             "expt_pka",
        "Expt. Uncertainty":     "expt_uncertainty",
        "Expt. Method":          "expt_method",
        "Reference":             "reference",
    })
    df["source_dataset"] = "PKAD2"
    df["uniprot_id"] = ""
    df["uniprot_resid"] = ""
    df["notes"] = ""
    # Normalise pdb_id to lowercase 4-char
    df["pdb_id"] = df["pdb_id"].str.strip().str.lower()
    df = df[df["pdb_id"].str.len() == 4].copy()
    df["chain"] = df["chain"].str.strip().str.upper()
    _remap = {"C-TERM": "CTE", "N-TERM": "NTE", "C-term": "CTE", "N-term": "NTE"}
    df["resname"] = df["resname"].str.strip().replace(_remap).str.upper()
    df["resid"] = df["resid"].str.strip()
    df["expt_pka"] = pd.to_numeric(df["expt_pka"], errors="coerce")
    df["expt_uncertainty"] = pd.to_numeric(df["expt_uncertainty"], errors="coerce")
    print(f"[PKAD2] {len(df)} rows after filtering")
    return df[_SCHEMA]


def read_pkad3_csv(path: Path) -> pd.DataFrame:
    """Read PKAD3.csv and return normalised DataFrame."""
    print(f"[PKAD3] Reading {path.name} …")
    df = pd.read_csv(path, dtype=str)
    # Drop duplicate header row (some PKAD3 exports repeat the header at row 2)
    df = df[df["id"] != "id"].copy()
    df = df.rename(columns={
        "pdb":          "pdb_id",
        "chain":        "chain",
        "resid":        "resid",
        "pKa":          "expt_pka",
        "uniprot_id":   "uniprot_id",
        "uniprot_resid":"uniprot_resid",
        "Ref":          "reference",
        "Notes":        "notes",
    })
    df["source_dataset"] = "PKAD3"
    df["expt_uncertainty"] = ""
    df["expt_method"] = ""
    # Map single-letter type → 3-letter resname
    df["resname"] = df["type"].str.strip().map(AA1_TO_3).fillna(df.get("type", ""))
    df["pdb_id"] = df["pdb_id"].str.strip().str.lower()
    df = df[df["pdb_id"].str.len() == 4].copy()
    df["chain"] = df["chain"].str.strip().str.upper()
    df["resid"] = df["resid"].str.strip()
    df["expt_pka"] = pd.to_numeric(df["expt_pka"], errors="coerce")
    df["uniprot_id"] = df["uniprot_id"].fillna("").str.strip()
    df["uniprot_resid"] = df["uniprot_resid"].fillna("").str.strip()
    df["reference"] = df["reference"].fillna("").str.strip()
    df["notes"] = df["notes"].fillna("").str.strip()
    print(f"[PKAD3] {len(df)} rows after filtering")
    return df[_SCHEMA]


def read_pkadR_csv(path: Path) -> pd.DataFrame:
    """Read PKAD-R.csv and return normalised DataFrame."""
    print(f"[PKAD-R] Reading {path.name} …")
    df = pd.read_csv(path, dtype=str)
    df = df.rename(columns={
        "PDB":              "pdb_id",
        "Chain":            "chain",
        "ResID in PDB":     "resid",
        "ResName":          "resname",
        "Expt. pKa":        "expt_pka",
        "Expt. Uncertainty":"expt_uncertainty",
        "Expt. Method":     "expt_method",
        "Reference":        "reference",
        "Notes":            "notes",
    })
    df["source_dataset"] = "PKAD-R"
    df["uniprot_id"] = ""
    df["uniprot_resid"] = ""
    # Normalise C-term / N-term labels
    remap = {"C-term": "CTE", "N-term": "NTE", "C-Term": "CTE", "N-Term": "NTE"}
    df["resname"] = df["resname"].str.strip().replace(remap).str.upper()
    df["pdb_id"] = df["pdb_id"].str.strip().str.lower()
    df = df[df["pdb_id"].str.len() == 4].copy()
    df["chain"] = df["chain"].str.strip().str.upper()
    df["resid"] = df["resid"].str.strip()
    df["expt_pka"] = pd.to_numeric(df["expt_pka"], errors="coerce")
    df["expt_uncertainty"] = pd.to_numeric(df["expt_uncertainty"], errors="coerce")
    df["reference"] = df["reference"].fillna("").str.strip()
    df["notes"] = df["notes"].fillna("").str.strip()
    print(f"[PKAD-R] {len(df)} rows after filtering")
    return df[_SCHEMA]


# ---------------------------------------------------------------------------
# 2. Super-table
# ---------------------------------------------------------------------------

def build_super_table(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Deduplicate across all three datasets by (pdb_id, chain, resname, resid).

    Adds boolean columns in_pkad2, in_pkad3, in_pkadR.
    When the same residue appears in multiple datasets with conflicting pKa values,
    stores all three in separate columns and sets pka_conflict=True.
    The merged expt_pka prefers PKAD3 > PKAD-R > PKAD2.
    """
    combined = pd.concat(frames, ignore_index=True)
    key = ["pdb_id", "chain", "resname", "resid"]

    # Per-dataset pKa columns
    ds_order = {"PKAD3": 0, "PKAD-R": 1, "PKAD2": 2}
    combined["_priority"] = combined["source_dataset"].map(ds_order).fillna(99)
    combined_sorted = combined.sort_values("_priority")

    rows = []
    for grp_key, grp in combined_sorted.groupby(key, sort=False):
        pdb_id, chain, resname, resid = grp_key
        datasets = sorted(grp["source_dataset"].unique(),
                          key=lambda d: ds_order.get(d, 99))

        pka_by_ds: dict[str, Optional[float]] = {}
        for _, r in grp.iterrows():
            ds = r["source_dataset"]
            v = r["expt_pka"]
            if pd.notna(v) and ds not in pka_by_ds:
                pka_by_ds[ds] = float(v)

        # Preferred merged pKa (first dataset in priority order that has a value)
        merged_pka = next(
            (pka_by_ds[d] for d in ("PKAD3", "PKAD-R", "PKAD2") if d in pka_by_ds),
            None,
        )

        pka_values = [v for v in pka_by_ds.values() if v is not None]
        pka_conflict = (len(pka_values) > 1 and
                        (max(pka_values) - min(pka_values)) > 0.1)

        # Take the first available value for non-pKa fields
        first = grp.iloc[0]
        uniprot_id = next(
            (r["uniprot_id"] for _, r in grp.iterrows()
             if r["uniprot_id"] and str(r["uniprot_id"]) not in ("", "nan")),
            "",
        )

        rows.append({
            "pdb_id":            pdb_id,
            "chain":             chain,
            "resname":           resname,
            "resid":             resid,
            "expt_pka":          merged_pka,
            "expt_pka_pkad2":    pka_by_ds.get("PKAD2"),
            "expt_pka_pkad3":    pka_by_ds.get("PKAD3"),
            "expt_pka_pkadR":    pka_by_ds.get("PKAD-R"),
            "pka_conflict":      pka_conflict,
            "expt_uncertainty":  first["expt_uncertainty"] if pd.notna(first["expt_uncertainty"]) else "",
            "uniprot_id":        uniprot_id,
            "uniprot_resid":     first["uniprot_resid"],
            "expt_method":       first["expt_method"],
            "reference":         first["reference"],
            "notes":             first["notes"],
            "in_pkad2":          "PKAD2" in datasets,
            "in_pkad3":          "PKAD3" in datasets,
            "in_pkadR":          "PKAD-R" in datasets,
        })

    df = pd.DataFrame(rows)
    df = df.sort_values(["pdb_id", "chain", "resid"]).reset_index(drop=True)
    print(f"[super_table] {len(df)} unique (pdb_id, chain, resname, resid) entries "
          f"from {len(combined)} raw rows")
    conflicts = df["pka_conflict"].sum()
    if conflicts:
        print(f"[super_table] {conflicts} entries have conflicting pKa values across datasets")
    return df


# ---------------------------------------------------------------------------
# 3. Crystal PDB download via Biopython
# ---------------------------------------------------------------------------

def download_crystal_pdbs(
    pdb_ids: list[str],
    out_dir: Path,
    skip_existing: bool = True,
) -> dict[str, bool]:
    """Download crystal PDB structures via Biopython PDBList, strip to heavy
    protein atoms, save to out_dir/<pdb_id>.pdb.

    Returns {pdb_id: success}.
    """
    PDBList, PDBParser, PDBIO, Select = _import_biopython()

    class _Sel(Select):
        def accept_model(self, m):
            return m.id == 0
        def accept_residue(self, r):
            return r.id[0] == " "  # standard protein residues only
        def accept_atom(self, a):
            elem = (a.element or "").strip().upper()
            if elem in ("H", "D"):
                return False
            if not elem and a.name.strip().startswith("H"):
                return False
            return True

    pdbl = PDBList(server="https://files.rcsb.org", verbose=False)
    parser = PDBParser(QUIET=True)
    results: dict[str, bool] = {}
    total = len(pdb_ids)

    for i, pdb_id in enumerate(sorted(pdb_ids), 1):
        out_file = out_dir / f"{pdb_id}.pdb"
        if skip_existing and out_file.exists():
            results[pdb_id] = True
            continue

        print(f"  [{i}/{total}] Downloading crystal {pdb_id} …", end=" ", flush=True)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                path = pdbl.retrieve_pdb_file(
                    pdb_id.upper(), file_format="pdb", pdir=tmp, overwrite=True
                )
                if not path or not Path(path).exists():
                    print("FAILED (no file)")
                    results[pdb_id] = False
                    continue
                structure = parser.get_structure(pdb_id, path)

            io = PDBIO()
            io.set_structure(structure)
            io.save(str(out_file), _Sel())
            print("ok")
            results[pdb_id] = True
        except Exception as exc:
            print(f"FAILED ({exc})")
            results[pdb_id] = False

    ok = sum(v for v in results.values())
    fail = len(results) - ok
    skip = total - len(results)
    print(f"[crystal PDBs] downloaded={ok}  failed={fail}  skipped(existing)={skip}")
    return results


# ---------------------------------------------------------------------------
# 4. AlphaFold DB download
# ---------------------------------------------------------------------------

def fetch_af3_pdb(
    uniprot_id: str,
    out_dir: Path,
    skip_existing: bool = True,
) -> tuple[bool, str]:
    """Download AF-{uniprot_id}-F1-model_v<N>.pdb from the AlphaFold DB.

    Tries v4 first, then v3, then v2.  Strips to heavy protein atoms.
    Returns (success, version_string).  Writes to out_dir/AF-{uniprot_id}-F1-model_v4.pdb
    (always uses the v4 filename regardless of which version was actually found).
    """
    canonical_out = out_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
    if skip_existing and canonical_out.exists():
        return True, "cached"

    for version in (4, 3, 2):
        url = (f"https://alphafold.ebi.ac.uk/files/"
               f"AF-{uniprot_id}-F1-model_v{version}.pdb")
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
            lines = _strip_heavy_protein_text(raw)
            canonical_out.write_text("\n".join(lines) + "\n")
            return True, f"v{version}"
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                continue
            raise
        except Exception as exc:
            raise RuntimeError(f"AF download failed for {uniprot_id}: {exc}") from exc

    return False, ""


def _strip_heavy_protein_text(pdb_text: str) -> list[str]:
    """Keep only ATOM lines with non-hydrogen elements."""
    out = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM  "):
            continue
        name = line[12:16].strip()
        elem = line[76:78].strip().upper() if len(line) > 76 else ""
        if elem in ("H", "D"):
            continue
        if not elem and name.startswith("H"):
            continue
        out.append(line)
    return out


def download_af3_pdbs(
    uniprot_pdb_pairs: list[tuple[str, str]],
    out_dir: Path,
    skip_existing: bool = True,
) -> tuple[dict[str, bool], list[dict]]:
    """Download AF3 predictions for all unique UniProt IDs.

    Parameters
    ----------
    uniprot_pdb_pairs: list of (uniprot_id, pdb_id) — used to label missing entries.
    out_dir: directory for AF3 PDB files.
    skip_existing: skip files already present.

    Returns
    -------
    results: {uniprot_id: success}
    missing: list of dicts {uniprot_id, pdb_id} for failed downloads.
    """
    seen: set[str] = set()
    results: dict[str, bool] = {}
    missing: list[dict] = []
    total = len({u for u, _ in uniprot_pdb_pairs})
    done = 0

    for uniprot_id, pdb_id in uniprot_pdb_pairs:
        if not uniprot_id or uniprot_id in seen:
            continue
        seen.add(uniprot_id)
        done += 1
        print(f"  [{done}/{total}] AF3 {uniprot_id} ({pdb_id}) …", end=" ")
        try:
            ok, ver = fetch_af3_pdb(uniprot_id, out_dir, skip_existing)
        except Exception as exc:
            print(f"ERROR ({exc})")
            ok, ver = False, ""
        results[uniprot_id] = ok
        if ok:
            print(f"ok ({ver})")
        else:
            print("not found")
            missing.append({"uniprot_id": uniprot_id, "pdb_id": pdb_id})

    ok_count = sum(v for v in results.values())
    print(f"[AF3 PDBs] downloaded={ok_count}  missing={len(missing)}")
    return results, missing


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdbs-only", action="store_true",
                        help="Skip table processing; only (re-)download PDB structures.")
    parser.add_argument("--no-af3", action="store_true",
                        help="Skip AlphaFold DB downloads.")
    parser.add_argument("--skip-existing", action="store_true", default=True)
    parser.add_argument("--no-skip-existing", dest="skip_existing", action="store_false")
    parser.add_argument("--data-dir", type=Path, default=HERE / "pkad_data",
                        help="Output data directory (default: pkad_data/).")
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    clean_dir = data_dir / "clean"
    crystal_dir = data_dir / "pdbs" / "heavy_atom"
    af3_dir = data_dir / "af3_pdbs"

    for d in (clean_dir, crystal_dir, af3_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: Read + normalise tables
    # ------------------------------------------------------------------
    if not args.pdbs_only:
        frames = []

        if PKAD2_SRC.exists():
            df2 = read_pkad2_ods(PKAD2_SRC)
            df2.to_csv(clean_dir / "PKAD2_clean.csv", index=False)
            frames.append(df2)
        else:
            print(f"[PKAD2] WARNING: source not found at {PKAD2_SRC}", file=sys.stderr)

        if PKAD3_SRC.exists():
            df3 = read_pkad3_csv(PKAD3_SRC)
            df3.to_csv(clean_dir / "PKAD3_clean.csv", index=False)
            frames.append(df3)
        else:
            print(f"[PKAD3] WARNING: source not found at {PKAD3_SRC}", file=sys.stderr)

        if PKADР_SRC.exists():
            dfR = read_pkadR_csv(PKADР_SRC)
            dfR.to_csv(clean_dir / "PKAD-R_clean.csv", index=False)
            frames.append(dfR)
        else:
            print(f"[PKAD-R] WARNING: source not found at {PKADР_SRC}", file=sys.stderr)

        if not frames:
            print("[ERROR] No PKAD tables found. Check repo root.", file=sys.stderr)
            sys.exit(1)

        super_table = build_super_table(frames)
        super_table.to_csv(data_dir / "super_table.csv", index=False)
        print(f"\n[write] {data_dir / 'super_table.csv'}")

    else:
        # Load existing super_table for PDB ID extraction
        super_table_path = data_dir / "super_table.csv"
        if not super_table_path.exists():
            print("[ERROR] super_table.csv not found. Run without --pdbs-only first.",
                  file=sys.stderr)
            sys.exit(1)
        super_table = pd.read_csv(super_table_path, dtype=str)

    # ------------------------------------------------------------------
    # Step 2: Crystal PDB downloads via Biopython
    # ------------------------------------------------------------------
    pdb_ids = sorted(super_table["pdb_id"].dropna().unique())
    print(f"\nDownloading {len(pdb_ids)} unique crystal PDB structures …")
    download_crystal_pdbs(pdb_ids, crystal_dir, skip_existing=args.skip_existing)

    # ------------------------------------------------------------------
    # Step 3: AlphaFold DB downloads (UniProt IDs from PKAD3)
    # ------------------------------------------------------------------
    if not args.no_af3:
        uniprot_col = super_table["uniprot_id"].fillna("")
        pkad3_mask = super_table.get("in_pkad3", pd.Series(False, index=super_table.index))
        if isinstance(pkad3_mask, pd.Series):
            pkad3_mask = pkad3_mask.astype(str).isin(["True", "true", "1"])
        has_uniprot = uniprot_col.str.strip().str.len() > 0

        af_rows = super_table[has_uniprot][["uniprot_id", "pdb_id"]].drop_duplicates(
            subset="uniprot_id"
        )
        pairs = list(zip(af_rows["uniprot_id"], af_rows["pdb_id"]))

        if pairs:
            print(f"\nDownloading {len(pairs)} unique AlphaFold DB predictions …")
            _, missing = download_af3_pdbs(pairs, af3_dir, skip_existing=args.skip_existing)
            # Write missing table
            missing_df = pd.DataFrame(missing, columns=["uniprot_id", "pdb_id"])
            missing_df.to_csv(data_dir / "missing_af3.csv", index=False)
            print(f"[write] {data_dir / 'missing_af3.csv'}  ({len(missing_df)} missing)")
        else:
            print("[AF3] No UniProt IDs found (PKAD3 not loaded?); skipping AF3 downloads.")

    print("\nDone.")


if __name__ == "__main__":
    main()
