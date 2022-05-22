redo-ifchange genes.label.do scripts/parse_seqs.py raw_data/dna_seqs.txt.gz
scripts/parse_seqs.py raw_data/dna_seqs.txt.gz -o genes
touch "$3"
