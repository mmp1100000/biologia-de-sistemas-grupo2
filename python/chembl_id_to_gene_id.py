import pickle

from chembl_webresource_client.new_client import new_client
from tqdm import tqdm
from python.drug_gene import get_tsv_gene_ids
from python.protein_entrez import get_gene_id
import fileinput


# chembl_file = '../outputs/interactions.tsv'


def get_accession_from_target(target):
    return target[0]['target_components'][0]['accession']


def chembl_id_to_accession_id():
    """
    Returns a dictionary of chembl_id:accession_id pairs, given the name of the interactions file
    @TODO FIND OUT WHY SOME ACCESSION IDS CAN'T BE FOUND
    :param interactions_file: interactions file
    :return: dict
    """
    target = new_client.target
    accession_ids = list()
    for chembl_id in tqdm(fileinput.input()):
        # chembl_id = line.split('\t')[0]
        try:
            res2 = target.filter(target_chembl_id=chembl_id)
            accession_ids.append(get_accession_from_target(res2))
        except Exception:
            accession_ids.append(" ")
    return accession_ids


accession_list = chembl_id_to_accession_id()

gene_list = get_gene_id(','.join(accession_list))  # Join the accession list into a single string

for gene_id in gene_list:
    print(gene_id)

#tsv_genes = get_tsv_gene_ids()  # Get gene ids from tsv

#for found_gene in gene_list:  # Compare found genes with tsv genes
#    if found_gene in tsv_genes:
#        print(found_gene)
