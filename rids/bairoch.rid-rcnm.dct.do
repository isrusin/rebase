redo-ifchange bairoch.rid-rcnm.dct.do ../raw_data/bairoch.txt.gz
zgrep -E '^(ID|AC)' ../raw_data/bairoch.txt.gz | tr -s '; ' '\t' | cut -f 2 |
    paste -s -d '\t\n' | sed -E 's/(.+)\tRB0*(.+)/\2\t\1/' | sort -u |
    sed '1i#:REBASE_ID\tComplex_name' > "$3"
