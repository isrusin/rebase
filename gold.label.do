redo-ifchange gold.label.do scripts/parse_seqs.py \
    raw_data/protein_gold_seqs.txt.gz
scripts/parse_seqs.py raw_data/protein_gold_seqs.txt.gz -o gold
touch "$3"
