redo-ifchange gold.label.do scripts/parse_seqs.py \
    raw_data/protein_gold_seqs.txt.gz
scripts/parse_seqs.py raw_data/protein_gold_seqs.txt.gz -o gold
gzip -f --rsyncable gold.fasta
touch "$3"
