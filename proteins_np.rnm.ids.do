redo-ifchange proteins_np.rnm.ids.do proteins.label
sieve_tsv.py -f7 -v 'no' proteins.tsv | cut -f1 > "$3"
