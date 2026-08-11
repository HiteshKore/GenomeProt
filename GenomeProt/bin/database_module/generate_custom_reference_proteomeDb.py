#!/usr/bin/env python3
"""
External helper script (NOT part of the core GenomeProt tool).
Use this script when you need to update a custom reference proteome
database in GenomeProt.

WHAT IT DOES
------------
Takes a UniProt and an OpenProt protein sequence database (both in
tabular format: "sequence_header<TAB>protein_sequence") and generates
a single, non-redundant combined database that includes following columns:
sequence	openprot_acc	accession	description	uniprot_reviewed_status	aggregated_accessions.

INPUT FILE SOURCES
-------------------
- OpenProt database:
    Download the database for your organism of interest from
    https://openprot.org/downloads

- UniProt database:
    Download the proteome that matches the OpenProt taxonomy ID
    ("TX=" field in the OpenProt header) from
    https://www.uniprot.org/proteomes

  Following UniProt proteome IDs were used generate custom proteome databases:
      UP000005640  Human
      UP000002277  Chimp
      UP000009136  Cow
      UP000000589  Mouse
      UP000002494  Rat
      UP000000803  Fruit fly
      UP000008143  Frog
      UP000000437  Zebrafish
      UP000002311  Yeast
      UP000001940  C. elegans

USAGE
-----
    python3 OpenProt_db_update.py <uniprot_db_tbl.fasta> <openprot_tbl.fasta> <ORGANISM_CODE> <outfile>

    Arguments:
        1. uniprot_db_tbl.fasta  - UniProt database in tabular format
        2. openprot_tbl.fasta    - OpenProt database in tabular format
        3. ORGANISM_CODE         - One of: HUMAN, MOUSE, CAEEL, DROME, RAT, DANRE
                                    (must match the organism token used in the
                                    UniProt FASTA header, e.g. ">sp|P12345|ACTB_HUMAN ...")
        4. outfile                - Path to write the combined, tab-separated output file

"""

import sys
import re
import pandas as pd

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
# args[1] = UniProt tabular database file
# args[2] = OpenProt tabular database file
# args[3] = Organism code (e.g. HUMAN, MOUSE, CAEEL, DROME, RAT, DANRE)
# args[4] = Output file path
args = sys.argv


# ---------------------------------------------------------------------------
# Step 1: Load the UniProt database
# ---------------------------------------------------------------------------
# Expected format: two tab-separated columns -> header, sequence (no header row)
uniprotdb = pd.read_csv(args[1], sep='\t', header=None)
uniprotdb.columns = ['header', 'sequence']


def parse_header(header, org):
    """
    Extract the gene/protein description from a UniProt FASTA header.

    UniProt headers typically look like:
        >sp|P12345|ACTB_HUMAN Actin, cytoplasmic 1 OS=Homo sapiens ...

    This function pulls out the description between
    the organism code (e.g. "HUMAN") and the "OS=" (organism species) tag.

    Parameters
    ----------
    header : str
        The full UniProt FASTA header.
    org : str
        The organism code.

    Returns
    -------
    str or None
        The extracted description, or None if the pattern is not found.
    """
    pattern = rf'{re.escape(org)} (.*?) OS='
    match = re.search(pattern, header)
    return match.group(1) if match else None


# Extract the UniProt accession from each header.
uniprotdb["accession"] = uniprotdb["header"].str.extract(r"^>\w+\|([^|]+)\|")

# Organism code provided on the command line (e.g. HUMAN, MOUSE, ...)
organism = args[3]

# Gene/protein description, parsed from the header using the organism code
uniprotdb['description'] = uniprotdb["header"].apply(parse_header, args=(organism,))

# UniProt reviewed/unreviewed status
# ">sp|" (reviewed) or ">tr|" (unreviewed).
uniprotdb['reviewed_status'] = uniprotdb['header'].str.split('|').str[0].str.lstrip('>')

# The raw header is not required once accession/description/status
# have been extracted from it.
uniprotdb = uniprotdb.drop(columns=["header"])

# ---------------------------------------------------------------------------
# Step 2: Collapse the UniProt table to unique sequences
# ---------------------------------------------------------------------------
# Multiple headers can map to the same protein sequence. Group by sequence
# and combine the associated accessions / descriptions / reviewed statuses
# into comma-separated, de-duplicated strings.
uniprotdb = (
    uniprotdb.groupby("sequence", as_index=False)
    .agg(
        accession=("accession", lambda x: ",".join(set(x))),
        description=("description", lambda x: ",".join(set(x))),
        uniprot_reviewed_status=("reviewed_status", lambda x: ",".join(sorted(set(x))))
    )
)


# ---------------------------------------------------------------------------
# Step 3: Load the OpenProt database
# ---------------------------------------------------------------------------
# Expected format: two tab-separated columns -> header, sequence
openprotdb = pd.read_csv(args[2], sep='\t', header=None)
openprotdb.columns = ['header', 'sequence']

# Extract the OpenProt accession from each header.
openprotdb["openprot_acc"] = openprotdb["header"].str.extract(r"^>([^|]+)")

# The raw header is not required once the accession has been extracted.
openprotdb = openprotdb.drop(columns=["header"])


# ---------------------------------------------------------------------------
# Step 4: Combine OpenProt and UniProt databases
# ---------------------------------------------------------------------------
# Full outer join on protein sequence: keeps sequences unique to OpenProt,
# unique to UniProt, and those found in both databases.
openprot_uniprot_combined = pd.merge(openprotdb, uniprotdb, on="sequence", how="outer")

# Replace missing/empty values with "-" as a placeholder for readability.
openprot_uniprot_combined["openprot_acc"] = openprot_uniprot_combined["openprot_acc"].fillna("-").replace("", "-")
openprot_uniprot_combined["accession"] = openprot_uniprot_combined["accession"].fillna("-").replace("", "-")
openprot_uniprot_combined["description"] = openprot_uniprot_combined["description"].fillna("-").replace("", "-")

# Build a single, de-duplicated, ordered list of accessions from both
# OpenProt and UniProt for each sequence.
openprot_uniprot_combined["aggregated_accessions"] = openprot_uniprot_combined.apply(
    lambda row: ",".join(
        dict.fromkeys(
            [x for x in (str(row["openprot_acc"]).split(",") +
                         str(row["accession"]).split(","))
             if x not in ["", "-", "nan", "NaN"]]
        )
    ),
    axis=1
)


# ---------------------------------------------------------------------------
# Step 5: Write the combined database to the output file
# ---------------------------------------------------------------------------
openprot_uniprot_combined.to_csv(args[4], sep='\t', index=False, header=True)
