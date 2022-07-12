#!/usr/bin/env bash

redo-ifchange default.nrnm.ids.do

if [[ "$2" =~ ^rebnr$ ]]; then
    redo-ifchange rebnr.clusters.tsv
    cut -f3 rebnr.clusters.tsv | sed '1s/^/#:/' | LC_ALL=C sort -u > "$3"
elif [[ "$2" =~ ^rebnr_np$ ]]; then
    redo-ifchange rebnr.nrnm.ids proteins_np.rnm.ids
    # If a REBNR cluster contains non-putatives, its representative will be
    # one of them.
    lists2list.py -i rebnr.nrnm.ids proteins_np.rnm.ids > "$3"
elif [[ "$2" =~ ^rebnr_gs$ ]]; then
    redo-ifchange proteins_gs.rnm.ids rebnr.clusters.tsv
    # If a REBNR cluster contains proteins from the gold standard,
    # the representative can be neither of them. But the cluster should
    # still be considered gold standard.
    cut -f 1,3 rebnr.clusters.tsv | replace_ids.py -ts - proteins_gs.rnm.ids |
        sed '1s/REBASE/REBNR/' | LC_ALL=C sort -u > "$3"
else
    echo "Bad target name $1!" >&2
    exit 1
fi
