#! /usr/bin/env python3

"""Parse protein names and use this info to fix their annotations.

REBASE contains some unconventional systems and proteins, but such
proteins are annotated as traditional R-M system components, except for
naming pattern. This script fixes type annotation of such proteins.
After the correction, it is possible (to some extent) to split proteins
by systems and complexes.
"""


import argparse
import logging
import re
import string
import sys
import gzip

from etsv import InputField, OutputField, ETSVReader, ETSVWriter


TYPE_RE = re.compile(
    """
    .*
    [^.-] # type can't be at the beginning of a name
    (?P<type>
        Mrr | Mcr[ABC] | Dam | Dcm | Gmr(?: S|D|SD) |
        Dnmt | CMT | DRM | MET | Dnd[ABCDE] |
        Dpt[FGH] | Ssp[ABCDFGH]
    )
    [0-9]* [A-Za-z]? (?: -.+)? P? $ # only some endings are allowed
    """, re.X
)

SYSTEM_TAGS = {
    "McrA": "McrABC", "McrB": "McrABC", "McrC": "McrABC",
    "GmrS": "GmrSD", "GmrD": "GmrSD",
    "DndA": "Dnd/Dpt", "DndB": "Dnd/Dpt", "DndC": "Dnd/Dpt",
    "DndD": "Dnd/Dpt", "DndE": "Dnd/Dpt",
    "DptF": "Dnd/Dpt", "DptG": "Dnd/Dpt", "DptH": "Dnd/Dpt",
    "SspA": "Ssp", "SspB": "Ssp", "SspC": "Ssp", "SspD": "Ssp",
    "SspF": "Ssp", "SspG": "Ssp", "SspH": "Ssp",
}

SYSTEM_TYPES = {
    "DndA": "Dnd/Dpt", "DndB": "Dnd/Dpt", "DndC": "Dnd/Dpt",
    "DndD": "Dnd/Dpt", "DndE": "Dnd/Dpt",
    "DptF": "Dnd/Dpt", "DptG": "Dnd/Dpt", "DptH": "Dnd/Dpt",
    "SspA": "Ssp", "SspB": "Ssp", "SspC": "Ssp", "SspD": "Ssp",
    "SspF": "Ssp", "SspG": "Ssp", "SspH": "Ssp",
    "M.Dcm": "Dcm", "V.Dcm": "Dcm",
    "Dam": "Orphan_M", "Dnmt": "Orphan_M", "CMT": "Orphan_M",
    "DRM": "Orphan_M", "MET": "Orphan_M",
}

CATEGORIES = {
    "Dnmt": "eM", "CMT": "eM", "DRM": "eM", "MET": "eM",
    "Mrr": "R", "McrA": "R", "McrB": "R", "McrC": "R",
    "V": "N", "V.Dcm": "N", "Nt": "N", "Nb": "N",
    "Dam": "M", "M.Dcm": "M",
    "GmrS": "GH", "GmrD": "GH", "GmrSD": "GH",
    "DndA": "PT", "DndB": "PT", "DndC": "PT", "DndD": "PT", "DndE": "PT",
    "DptF": "PT", "DptG": "PT", "DptH": "PT",
    "SspA": "PT", "SspB": "PT", "SspC": "PT", "SspD": "PT",
    "SspF": "PT", "SspG": "PT", "SspH": "PT",
}


def add_type_tag(entry):
    name = entry["name"]
    if "." in name:
        prot_type = entry["prot_type"]
        type_tag = name.split(".", 1)[0].rstrip("0123456789")
        if type_tag != prot_type:
            logging.info(
                f"{name} tag does not match protein type {prot_type}"
            )
        entry["type_tag"] = type_tag


def add_name_type(entry):
    type_m = TYPE_RE.match(entry["name"])
    if type_m:
        name_type = type_m.group("type")
        if name_type == "Dcm":
            name_type = f"{entry['prot_type']}.Dcm"
        entry["name_type"] = name_type


def fix_ssps(entries):
    ssp_types = dict()
    ssp_prots = dict()
    for entry in entries:
        if entry["sys_type"] != "Type II":
            logging.info(
                f"{entry['name']} is not Ssp, system type mismatch"
            )
            del entry["name_type"]
        elif entry["prot_type"] not in "MR":
            logging.info(
                f"{entry['name']} is not Ssp, protein type mismatch"
            )
            del entry["name_type"]
        else:
            name = entry["name"]
            sys, type_tag = name.rsplit("Ssp", 1)
            ssp_types.setdefault(sys, set()).add(type_tag[0])
            ssp_prots.setdefault(sys, []).append(entry)
    for sys, types in ssp_types.items():
        if len(types) == 1:
            for entry in ssp_prots[sys]:
                logging.info(
                    f"{entry['name']} is not Ssp, no other components"
                )
                del entry["name_type"]


def add_system(entry):
    name = entry["name"]
    complex_name = entry["complex"]
    if complex_name:
        name = complex_name
    if name.endswith("P"):
        name = name[:-1]
    if "type_tag" in entry:
        _tag, name = name.split(".", 1)
    name_type = entry.get("name_type")
    if name_type in SYSTEM_TAGS:
        name = SYSTEM_TAGS[name_type].join(name.rsplit(name_type, 1))
    entry["system"] = name


def add_category(entry):
    prot_type = entry["prot_type_fixed"]
    entry["category"] = CATEGORIES.get(prot_type, prot_type)


def fix_system_type(entry):
    name_type = entry.get("name_type", "None")
    if name_type in SYSTEM_TYPES:
        entry["sys_type_fixed"] = SYSTEM_TYPES[name_type]
        logging.info(f"change {entry['name']} system type")
    else:
        entry["sys_type_fixed"] = entry["sys_type"].replace(" ", "_")


def fill_and_write(entry, outsv):
    if "name_type" in entry:
        entry["prot_type_fixed"] = entry["name_type"]
    elif "type_tag" in entry:
        entry["prot_type_fixed"] = entry["type_tag"]
    else:
        entry["prot_type_fixed"] = entry["prot_type"]
    add_category(entry)
    add_system(entry)
    fix_system_type(entry)
    outsv.write_entry(entry)


def assign_systems(entries, outsv):
    ssps = []
    # There are errors in REBASE with Ssp component names,
    # store them for futher analysis

    for entry in entries:
        add_type_tag(entry)
        add_name_type(entry)
        if "Ssp" in entry.get("name_type", "None"):
            ssps.append(entry)
            logging.info(
                f"keep {entry['name']} for futher Ssp analysis"
            )
        else:
            fill_and_write(entry, outsv)

    fix_ssps(ssps)
    for entry in ssps:
        fill_and_write(entry, outsv)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Parse REBASE names.")
    parser.add_argument(
        "intsv", metavar="TSV", type=argparse.FileType("r"),
        help="input TSV file with REBASE proteins"
    )
    parser.add_argument(
        "-o", "--outsv", metavar="TSV", type=argparse.FileType("w"),
        default=sys.stdout, help="output TSV, default is STDOUT"
    )
    args = parser.parse_args(argv)
    input_fields = [
        InputField("name", "REBASE_name"),
        InputField("complex", "Complex_name"),
        InputField("prot_type", "Protein_type"),
        InputField("sys_type", "System_type"),
    ]
    output_fields = [
        OutputField("name", "REBASE_name"),
        OutputField("complex", "Complex_name"),
        OutputField("prot_type", "Protein_type"),
        OutputField("prot_type_fixed", "Protein_type_fixed"),
        OutputField("category", "Protein_category"),
        OutputField("system", "System_name"),
        OutputField("sys_type", "System_type"),
        OutputField("sys_type_fixed", "System_type_fixed"),
    ]
    with args.intsv as intsv, args.outsv as outsv:
        assign_systems(
            ETSVReader(intsv, input_fields), ETSVWriter(outsv, output_fields)
        )


if __name__ == "__main__":
    logging.basicConfig(
        format="%(levelname)s: %(message)s", level=logging.INFO
    )
    sys.exit(main())

