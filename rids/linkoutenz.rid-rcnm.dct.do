redo-ifchange ../raw_data/linkoutenz.label
zcat ../raw_data/linkoutenz*.xml.gz |
    sed -nE 's/.*<Rule>([0-9]+).*/\1/p;s/.*<UrlName.*enzyme ([^<]+).*/\1/p' |
    paste -s -d '\t\n' | sort -u -k1,1n | sed '1i#:REBASE_ID\tComplex_name' \
    > "$3"
