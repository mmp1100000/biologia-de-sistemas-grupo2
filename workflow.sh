#! /usr/bin/env bash

python python/gsk_filter.py | sort -u > outputs/interactions.tsv

cut -f 1 < outputs/interactions.tsv | python python/chembl_id_to_gene_id.py