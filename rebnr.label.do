redo-ifchange rebnr.label.do scripts/select_uniq.py proteins.label
scripts/select_uniq.py proteins.fasta.gz \
    -o rebnr -n 'REBASE_name' -r 'REBNR_name'
touch "$3"
