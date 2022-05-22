redo-ifchange linkoutenz.label.do VERSION

for file in $(curl -sl "ftp.neb.com/pub/rebase/" | grep '^linkoutenz'); do
    curl -s "ftp://ftp.neb.com/pub/rebase/${file}" |
    gzip --rsyncable > "${file}.gz"
done

touch "$3"
