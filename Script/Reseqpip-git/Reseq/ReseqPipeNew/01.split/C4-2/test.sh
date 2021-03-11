for ID in `while read line; do echo $line|tr -d '\n'; done < $1`; do
    echo $ID
    echo '1'
done
