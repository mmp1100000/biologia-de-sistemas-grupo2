#! /usr/bin/env bash

cd `dirname $0`

prettyvis=false

# Parse options
while getopts ":v" opt; do
    case $opt in
        v)
            prettyvis=true
            ;;
        \?)
            echo "Invalid option: -$OPTARG."
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))

if [ $prettyvis = true ]; then
    (head -n 1 data/other_data_for_PKIS_compounds.csv | cut -d, -f 9,10,11,12,13,14,16,17 | tr ',' '\t' \
    && grep $1 < data/other_data_for_PKIS_compounds.csv | \
    awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '{print $9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$16"\t"$17}') | \
    sed 's/\t\t/\t \t/g;s/\t\t/\t \t/g' | column -s $'\t' -t | less
    exit 0
fi

head -n 1 data/other_data_for_PKIS_compounds.csv | cut -d, -f 9,10,11,12,13,14,16,17 | tr ',' '\t' \
&& grep $1 < data/other_data_for_PKIS_compounds.csv | \
awk -vFPAT='([^,]*)|("[^"]+")' -vOFS=, '{print $9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$16"\t"$17}'
