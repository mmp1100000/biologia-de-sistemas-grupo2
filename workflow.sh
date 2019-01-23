#! /usr/bin/env bash

set -e
cd `dirname $0`

printinfo()
{
    echo -e "\e[36m\e[1m$1\e[0m"
}

printerror()
{
    echo -e "\e[31m\e[1m$1\e[0m"
}

skip_gsk_filter=false
confidence_interval=0.95
min_inhibition_percentage=50

mkdir -p outputs

# Parse options
while getopts ":c:i:s" opt; do
    case $opt in
        s)
            skip_gsk_filter=true
            ;;
        c)
            confidence_interval=$OPTARG
            ;;
        i)
            min_inhibition_percentage=$OPTARG
            ;;
        \?)
            printerror "Invalid option: -$OPTARG."
            exit 1
            ;;
        :)
            printerror "Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))

if [ $skip_gsk_filter = false ]; then
    printinfo "Downloading unique interactions with minimum inhibition percentage of $min_inhibition_percentage%."
    python python/gsk_filter.py $min_inhibition_percentage | sort -u > outputs/interactions.tsv
fi

printinfo "Translating target CHEMBL ids to NCBI gene ids."
cut -f 1 < outputs/interactions.tsv | sort -u > outputs/unique_targets_chembl.tsv
python python/chembl_id_to_accession.py < outputs/unique_targets_chembl.tsv > outputs/unique_targets_accession.tsv
python python/accession_to_gene_id.py < outputs/unique_targets_accession.tsv > outputs/unique_targets_geneid.tsv
paste outputs/unique_targets_chembl.tsv outputs/unique_targets_geneid.tsv > outputs/targets_chembl_geneid.tsv
awk 'FNR==NR { a[$1]=$2; next } a[$1] { print $2"\t"a[$1] }' outputs/targets_chembl_geneid.tsv outputs/interactions.tsv > outputs/translated_interactions.tsv

printinfo "Parsing known breast cancer genes."
grep breast < data/journal.pcbi.1004120.s004.TSV | awk '{print $3}' | tr '/' '\n' > outputs/breast_cancer_genes.tsv

printinfo "Checking off breast cancer gene interactions."
awk 'FNR==NR { a[$0]; next } { if ($2 in a) { print $0"\t"1 } else { print $0"\t"0 } }' \
outputs/breast_cancer_genes.tsv outputs/translated_interactions.tsv | sort -k 3 -r > outputs/flagged_interactions.tsv

printinfo "Performing the overrepresentation statistical test."
python python/binom.py outputs/flagged_interactions.tsv $confidence_interval

printinfo "Generating graph."
python python/graph.py outputs/pvalued_interactions_${confidence_interval}.tsv $confidence_interval