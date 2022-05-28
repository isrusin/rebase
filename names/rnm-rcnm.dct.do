redo-ifchange rnm-rcnm.dct.do name_data.tsv
cut -f 1,2 name_data.tsv | sieve_tsv.py -rv '' -f2 > "$3"
