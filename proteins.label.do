redo-ifchange proteins.label.do scripts/parse_seqs.py \
    raw_data/protein_org_seqs.txt.gz
scripts/parse_seqs.py raw_data/protein_org_seqs.txt.gz -o proteins
gzip -f --rsyncable proteins.fasta
touch "$3"
