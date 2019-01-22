#! /usr/bin/env bash

# python python/gsk_filter.py | sort -u > outputs/interactions.tsv

cut -f 1 < outputs/interactions.tsv | sort -u > outputs/unique_targets_chembl.tsv

python python/chembl_id_to_accession.py < outputs/unique_targets_chembl.tsv > outputs/unique_targets_accession.tsv

python python/accession_to_gene_id.py < outputs/unique_targets_accession.tsv > outputs/unique_targets_geneid.tsv

paste outputs/unique_targets_chembl.tsv outputs/unique_targets_geneid.tsv > outputs/targets_chembl_geneid.tsv

awk 'FNR==NR{ a[$1]=$2; next} a[$1]{ print $2"\t"a[$1]}' outputs/targets_chembl_geneid.tsv outputs/interactions.tsv > outputs/translated_interactions.tsv
