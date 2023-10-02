redo-always
mkdir -p ../pages/
find ../pages/ -iname '*.html' | xargs grep -h -o 'TITLE.*TITLE' |
    sed 's/Enz /\t/;s/ - /\t/;s/<\//\t/' | cut -f 2,3 | sort -u |
    sed '1i#:REBASE_ID\tComalex_name' > "$3"
redo-stamp < "$3"
