IndexError analysis — 1hpx and 3fx5
Investigated 2026-06-23 from logs and TAPBS output on flashnas2.

Summary
=======

The "list index out of range" crash is a symptom of a TAPBS failure, not a
Python bug.  TAPBS terminates abnormally because of an unknown atom in the NTE
(N-terminal patch) state file, which leaves Karlsberg with no input to read,
which leaves the code with an empty result dict that it then tries to index.
Three things go wrong in sequence:

  1. TAPBS dies on HT1 in NTE (root cause, physics/data issue)
  2. Karlsberg is called anyway with a missing .pkint file (silent error)
  3. Code crashes on an empty occupancy dict (code robustness bug)


Step-by-step failure chain
==========================

Step 1 — TAPBS abnormal termination
------------------------------------
TAPBS reads all titratable residue state files (.st) to build the interaction
network.  For both 1hpx and 3fx5, it processes all residues successfully until
it reaches NTE-1_A:

    NTE-1_A: parsing NTE.st ... finished.
    Abnormal termination due to fatal error:
    Unknown atoms in NTE-1_A: R
    Unknown atom named HT1 by State::check. (Error 14)

NTE is the N-terminal patch (applied to residue 1 of each chain).  The NTE.st
state file lists atom HT1 as a protonatable proton.  TAPBS's internal atom-type
table does not recognise "HT1" (it reports the type as "R", an unknown code).

As a result, TAPBS exits without writing its normal output files:
  — init_calc.pkint    (pKint values per residue, input for Karlsberg)
  — init_calc.g        (pairwise interaction matrix, input for Karlsberg)
Only the intermediate .st parsing output and init_calc.sites are present.

Both structures are homodimers (chains A and B), which is why TAPBS processes
all the chain-B residues fine and only hits the error when it reaches NTE-1_A
at the end.  Chain-B NTE is absent because NTE is only applied once per chain
and chain B's NTE would be NTE-1_B.  Looking at the TAPBS output, it appears
NTE-1_B is not in the titratable list (only NTE-1_A is listed), so only one
NTE is encountered.

Step 2 — Karlsberg called with missing input
---------------------------------------------
The pKa calculation code does not check whether TAPBS succeeded before calling
Karlsberg.  Karlsberg is launched with a reference to the missing .pkint:

    conformation .../init_calc.pkint .../init_calc.g 0.00 kJ/mol

Karlsberg immediately terminates:

    I am sorry, abnormal termination with the following message:
    File '.../init_calc.pkint' given in line number 1 not found or
    not readable. (Error 6)

Karlsberg finishes in < 1 second with no occupancy output.

Step 3 — Empty dict IndexError
--------------------------------
parse_karlsberg_results() reads karlsberg.out and finds no occupancy data:

    Conformer occupancies: {}
    Occupancy keys: dict_keys([])

Downstream code then attempts to access the first element of this empty dict
(pattern: list(occupancies.values())[0] or equivalent), which raises:

    IndexError: list index out of range

This propagates up as "calc_pkas failed: list index out of range".


Root cause — HT1 in NTE.st
============================

The NTE patch (N-terminal) adds protonation states for the alpha-amine.  The
state file NTE.st lists HT1 as one of the atoms that changes between
protonation states.  HT1 is the C22/C36 atom name for the first amine proton
on the N-terminus.  TAPBS's own atom-type table appears not to include HT1 by
name; it resolves to type "R" (unknown).

This is likely a mismatch between the TAPBS .st file definition and the TAPBS
binary's internal atom table.  The NTE.st file probably needs to be updated to
use the atom name that TAPBS's binary expects.

Why only 1hpx and 3fx5?
  Most structures don't have a titratable NTE in the experimental dataset —
  it's only included when there's a measured pKa for the N-terminus.  Both 1hpx
  and 3fx5 have an experimental NTE pKa in the PKAD dataset, so the NTE
  residue is included in the TAPBS calculation for them.  No other structure in
  the benchmark has both a titratable NTE and completed TAPBS.

  1hpx: HIV-1 protease dimer — NTE-1_A titratable
  3fx5: another dimer       — NTE-1_A titratable


What needs to be fixed
=======================

Fix A (root cause) — update NTE.st
  Inspect the NTE.st state file used by TAPBS and find the correct atom name
  for HT1 according to the TAPBS binary's internal table.  This likely requires
  checking what atom names the TAPBS binary accepts for N-terminal amine protons
  (possibly "H" or "HN" instead of "HT1") and updating NTE.st accordingly.
  Location: wherever the .st files are stored in the KB3 installation.

Fix B (code robustness) — check Karlsberg output before indexing
  parse_karlsberg_results() should detect an empty occupancy dict and raise a
  clear, descriptive exception rather than returning an empty object that later
  crashes on first access.  A check like:
    if not occupancies:
        raise RuntimeError("Karlsberg returned no occupancy data — "
                           "check TAPBS output for fatal errors")
  would turn the opaque IndexError into an actionable error message.

Fix C (code robustness) — check TAPBS exit status
  The code should verify that TAPBS produced init_calc.pkint before proceeding
  to call Karlsberg.  If the file is absent, raise immediately with the last
  lines of init_calc.tapbs.out for context.


Resubmission note
==================

Both 1hpx and 3fx5 will fail again on resubmission until Fix A (NTE.st atom
name) is resolved.  They are not candidates for --ph7-only either, because
TAPBS runs for all three PACs and would hit the same error at the pH-7 PAC.
