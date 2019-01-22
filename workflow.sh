#! /usr/bin/env bash

python python/gsk_filter.py | sort -u > outputs/interactions.tsv

cut -f 1 < outputs/interactions.tsv | sort -u > outputs/unique_targets_chembl.tsv

python python/chembl_id_to_accession.py < outputs/unique_targets_chembl.tsv > outputs/unique_targets_accession.tsv

python python/accession_to_gene_id.py < outputs/unique_targets_accession.tsv > outputs/unique_targets_geneid.tsv

paste outputs/unique_targets_chembl.tsv outputs/unique_targets_geneid.tsv > outputs/targets_chembl_geneid.tsv

awk 'FNR==NR{ a[$1]=$2; next} a[$1]{ print $2"\t"a[$1]}' outputs/targets_chembl_geneid.tsv outputs/interactions.tsv > outputs/translated_interactions.tsv

grep breast < data/journal.pcbi.1004120.s004.TSV | awk '{print $3}' | tr '/' '\n' > outputs/breast_cancer_genes.tsv

awk 'FNR==NR{a[$0]; next} { if ($2 in a) { print $0"\t"1 } else { print $0"\t"0 } }' \
outputs/breast_cancer_genes.tsv outputs/translated_interactions.tsv | sort -k 3 -r > outputs/flagged_interactions.tsv
