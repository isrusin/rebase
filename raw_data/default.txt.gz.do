redo-ifchange default.txt.gz.do VERSION
curl -s "ftp://ftp.neb.com/pub/rebase/$2.txt" | gzip --rsyncable > "$3"
