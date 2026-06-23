"""
Strip hydrogen atoms from PDB files, writing heavy-atom-only output.

Handles:
  - ATOM / HETATM records
  - Hydrogen identification by element column (col 77-78) or atom name
  - Deuterium (D) treated same as hydrogen
  - Input can be a single PDB file or a directory of PDB files
  - Renumbers ATOM serial numbers after stripping

Usage:
    # Strip one file
    python strip_hydrogens.py input.pdb output.pdb

    # Strip all PDBs in a directory
    python strip_hydrogens.py pdbs/heavy_atom/ stripped/

    # Strip and overwrite in-place (creates backup .bak)
    python strip_hydrogens.py input.pdb --inplace

Note:
    The crystal structures in pkad_data/pdbs/heavy_atom/ are ALREADY stripped
    of hydrogens (done by download_pkad_data.py during data preparation).
    This script is provided for stripping new or user-supplied structures.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path


# Atom names that unambiguously indicate hydrogen even without element column
# KB3 and CHARMM hydrogen naming conventions
_H_NAME_PREFIXES = ("H", "D", "1H", "2H", "3H", "4H", "1D", "2D", "3D")


def _is_hydrogen(line: str) -> bool:
    """Return True if the ATOM/HETATM line is a hydrogen atom."""
    if len(line) < 14:
        return False

    # PDB element column: cols 77-78 (0-indexed: 76-77), right-justified
    if len(line) >= 78:
        element = line[76:78].strip().upper()
        if element in ("H", "D"):
            return True
        if element and element not in ("", " "):
            # Element is something else (not blank) — definitely not hydrogen
            return False

    # Fallback: atom name in cols 13-16 (0-indexed: 12-15)
    atom_name = line[12:16].strip().upper()
    # Strip leading digits (e.g. '1HB' -> 'HB')
    clean = atom_name.lstrip("0123456789")
    if clean.startswith(("H", "D")):
        return True

    return False


def strip_file(input_path: Path, output_path: Path, keep_ter: bool = True) -> tuple:
    """
    Strip hydrogens from a single PDB file.

    Returns (n_total, n_stripped, n_kept) atom record counts.
    """
    n_total = 0
    n_stripped = 0
    output_lines = []
    serial = 0

    with open(input_path) as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            record = line[:6] if len(line) >= 6 else line

            if record in ("ATOM  ", "HETATM"):
                n_total += 1
                if _is_hydrogen(line):
                    n_stripped += 1
                    continue
                # Renumber serial
                serial += 1
                # Rebuild serial field (cols 7-11, 1-indexed: 6-10)
                new_serial = f"{serial:5d}"
                line = line[:6] + new_serial + line[11:]
                output_lines.append(line + "\n")

            elif record == "TER   " and keep_ter:
                output_lines.append(raw_line)

            elif record not in ("ATOM  ", "HETATM", "TER   "):
                # Keep all non-atom records (HEADER, REMARK, CONECT, END, etc.)
                output_lines.append(raw_line)

    n_kept = n_total - n_stripped

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.writelines(output_lines)

    return n_total, n_stripped, n_kept


def strip_directory(input_dir: Path, output_dir: Path, verbose: bool = True) -> dict:
    """
    Strip all *.pdb files in input_dir, writing to output_dir.

    Returns a summary dict: pdb_id -> {n_total, n_stripped, n_kept}.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {}

    pdb_files = sorted(input_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"WARNING: no .pdb files found in {input_dir}", file=sys.stderr)
        return summary

    for pdb_file in pdb_files:
        out_file = output_dir / pdb_file.name
        n_total, n_stripped, n_kept = strip_file(pdb_file, out_file)
        pdb_id = pdb_file.stem
        summary[pdb_id] = {"n_total": n_total, "n_stripped": n_stripped, "n_kept": n_kept}
        if verbose:
            print(f"  {pdb_id}: {n_total} atoms → stripped {n_stripped} H → {n_kept} heavy atoms  →  {out_file}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", type=Path, help="Input PDB file or directory")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Output PDB file or directory (default: derived from input)",
    )
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite input file in-place (creates .bak backup first)",
    )
    parser.add_argument(
        "--no-ter",
        action="store_true",
        help="Drop TER records from output",
    )
    parser.add_argument("--quiet", "-q", action="store_true")
    args = parser.parse_args()

    keep_ter = not args.no_ter
    verbose = not args.quiet

    if args.input.is_dir():
        # Directory mode
        out_dir = args.output if args.output else args.input.parent / (args.input.name + "_heavy")
        if verbose:
            print(f"Stripping PDB directory: {args.input} → {out_dir}")
        summary = strip_directory(args.input, out_dir, verbose=verbose)
        if verbose:
            total_in = sum(v["n_total"] for v in summary.values())
            total_out = sum(v["n_kept"] for v in summary.values())
            print(f"\nTotal: {len(summary)} files, {total_in} input atoms → {total_out} heavy atoms")

    elif args.input.is_file():
        if args.inplace:
            backup = args.input.with_suffix(".pdb.bak")
            shutil.copy2(args.input, backup)
            if verbose:
                print(f"Backup: {backup}")
            out_path = args.input
        else:
            if args.output:
                out_path = args.output
            else:
                out_path = args.input.with_stem(args.input.stem + "_heavy")

        n_total, n_stripped, n_kept = strip_file(args.input, out_path, keep_ter=keep_ter)
        if verbose:
            print(f"{args.input.name}: {n_total} atoms → stripped {n_stripped} H → {n_kept} heavy atoms")
            print(f"Output: {out_path}")

    else:
        print(f"ERROR: input not found: {args.input}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
