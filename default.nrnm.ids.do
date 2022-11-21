#!/usr/bin/env bash

redo-ifchange default.nrnm.ids.do

if [[ "$2" =~ ^rebnr$ ]]; then
    redo-ifchange rebnr.clusters.tsv
    cut -f3 rebnr.clusters.tsv | sed '1s/^/#:/' | LC_ALL=C sort -u > "$3"
elif [[ "$2" =~ ^rebnr_(np|gs)$ ]]; then
    pset="${BASH_REMATCH[1]}"
    redo-ifchange rebnr.nrnm.ids "proteins_${pset}.rnm.ids"
    # If a REBNR cluster contains Gold Standard members or non-putatives,
    # its representative will be one of them, because of selection procedure.
    lists2list.py -i rebnr.nrnm.ids "proteins_${pset}.rnm.ids" > "$3"
else
    echo "Bad target name $1!" >&2
    exit 1
fi
