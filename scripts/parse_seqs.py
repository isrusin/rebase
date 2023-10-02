#! /usr/bin/env python3

"""Parse REBASE sequence files.

The REBASE sequence format is very similar to Fasta, except the following:
 1) each sequence ends with '<>' symbols;
 2) title lines have an internal structure of tab-separated key:value pairs.
The keys differ between different files, some of them are obligatory for all
entires, and others are not.

This script parses an input file in the REBASE sequence format and creates
separate fasta (only sequences with names as IDs) and TSV (all non-sequence
data) files.

Each description key corresponds to a column in the ouput TSV file.
Several keys have column positions determined (for some of them, values are
even further processed), all other columns arranged alphabetically by key
names.

Internally parsed keys include:
 * REBASE --> "REBASE_name", "Complex_name"
 * GenBank --> "Sequence_AC", "Sequence_source"
 * EnzType --> "System_type", "Protein_type", "Putative"
"""


import argparse
import logging
import re
import string
import sys
import gzip


TYPE_RE = re.compile(
    """
    (?: (?P<putative> putative) \s)?
    (?: (?P<systype> Type \s [IVG]+) \s)?
    (?: (?P<md> methyl-directed) \s)?
    (?P<activity>.+)
    """, re.X
)

NAME_RE = re.compile(
    """
    (?P<complex> [^ ]+)
    (?: \s \( (?P<name>[^ ]+) \))?
    (?: \s \( (?P<_altrm>RM\..+) \))?
    """, re.X
)

ACTIVITIES = {
    "methyltransferase": "M",
    "restriction enzyme": "R",
    "restriction enzyme/methyltransferase": "RM",
    "specificity subunit": "S",
    "control protein": "C",
    "nicking endonuclease": "V",
    "helicase domain protein": "H",
    "orphan methyltransferase": "M_orphan",
    "homing endonuclease": "R_homing",
    "methyltransferasespecificity subunit": "M",
}

COLUMNS = [
    "REBASE_name", "Complex_name",
    "Sequence_AC", "Sequence_source",
    "System_type", "Protein_type", "Putative",
]

RESIDUES = set(string.ascii_uppercase)


def parse_EnzType(type_str):
    type_match = TYPE_RE.match(type_str)
    assert type_match, f"Can't parse EnzType \"{type_str}\"!"

    is_put = "yes" if type_match.group("putative") else "no"
    sys_tp = type_match.group("systype") or "-"
    is_md = type_match.group("md")
    activity = type_match.group("activity").strip()
    assert activity in ACTIVITIES, f"Unknown enzyme type \"{activity}\"!"
    if "methyltransferasespecificity" in activity:
        logging.warning(f"strange activity '{activity}', treat as MTase")

    enz_tp = ACTIVITIES[activity]
    if enz_tp == "M_orphan":
        sys_tp, enz_tp = "Orphan M", "M"
    if enz_tp == "R_homing":
        sys_tp, enz_tp = "Homing", "R"
    if is_md and sys_tp == "Type II":
        sys_tp = "Type IIM"

    return {"Putative": is_put, "Protein_type": enz_tp, "System_type": sys_tp}


def parse_AC(ac_str):
    source = "INSDC"
    if ac_str.startswith("NEB"):
        source = "NEB"
    elif "_" in ac_str:
        source = "RefSeq"
    return {"Sequence_AC": ac_str, "Sequence_source": source}


def parse_Name(name_str):
    name_match = NAME_RE.match(name_str)
    assert name_match, f"Can't parse enzyme name {name_str}!"
    complex_name = name_match.group("complex")
    name = name_match.group("name")
    if name is None:
        name, complex_name = complex_name, ""
    elif name.startswith("RM.") or (abs(len(name) - len(complex_name)) > 2):
        logging.info(f"{complex_name} has alternative name {name}")
        name, complex_name = complex_name, ""
    return {"REBASE_name": name, "Complex_name": complex_name}


def write_seq(oufasta, name, seq_list):
    list_prepared = []
    for seq_frag in seq_list:
        frag_prepared = seq_frag.strip().replace("<>", "").replace(" ", "")
        if set(frag_prepared) - RESIDUES:
            logging.warning(f"skip {name} bad sequence line: '{seq_frag}'")
            continue
        if frag_prepared:
            list_prepared.append(frag_prepared)
    if list_prepared:
        print(">", name, sep="", file=oufasta)
        print(*list_prepared, sep="\n", file=oufasta)
    else:
        logging.warning(f"{name} has no sequence")


def write_entries(outsv, entries):
    column_set = set()
    for pairs in entries.values():
        column_set.update(pairs)
    columns_raw = sorted(column_set - set(COLUMNS))
    with outsv:
        title = COLUMNS + [f"{col}_raw" for col in columns_raw]
        print("#:", "\t".join(title), sep="", file=outsv)
        columns = COLUMNS + columns_raw
        for name, pairs in sorted(entries.items()):
            print(
                name, *(pairs.get(col, "") for col in columns[1:]),
                sep="\t", file=outsv
            )


def parse_seqs(inseqs, oufasta, outsv):
    entries = dict()
    with inseqs, oufasta:
        name = None
        seq_list = None
        for line in inseqs:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    write_seq(oufasta, name, seq_list)
                seq_list = []
                pairs = dict()
                for pair in line[1:].split("\t"):
                    key, value = map(str.strip, pair.split(":", 1))
                    pairs[key] = value
                pairs.update(parse_Name(pairs.pop("REBASE")))
                pairs.update(parse_EnzType(pairs.pop("EnzType")))
                pairs.update(parse_AC(pairs.pop("GenBank")))
                name = pairs.pop("REBASE_name")
                assert name not in entries, f"Repeated REBASE name {name}!"
                entries[name] = pairs
            elif name:
                seq_list.append(line)
        write_seq(oufasta, name, seq_list)
    write_entries(outsv, entries)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Parse REBASE seqs.txt file.")
    parser.add_argument(
        "inseqs", metavar="FILE",
        help="input REBASE *seqs.txt file, can be gzipped"
    )
    parser.add_argument(
        "-o", "--output-tag", dest="outag", metavar="TAG",
        help="""common TAG for the output files: TAG.fasta and TAG.tsv will
        be created; default TAG is the name of the input file"""
    )
    args = parser.parse_args(argv)
    inseqs_name = args.inseqs
    if inseqs_name.endswith(".gz"):
        inseqs = gzip.open(inseqs_name, mode="rt")
    else:
        inseqs = open(inseqs_name)
    outag = args.outag
    if outag is None:
        outag = inseqs.name
    parse_seqs(
        inseqs, open(f"{outag}.fasta", "w"), open(f"{outag}.tsv", "w")
    )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(levelname)s: %(message)s", level=logging.INFO
    )
    sys.exit(main())

