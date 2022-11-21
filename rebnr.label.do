redo-ifchange rebnr.label.do scripts/select_uniq.py \
    proteins_gs.rnm.ids proteins.label

scripts/select_uniq.py proteins.fasta.gz \
    -o rebnr -n 'REBASE_name' -r 'REBNR_name' --gold proteins_gs.rnm.ids
touch "$3"
