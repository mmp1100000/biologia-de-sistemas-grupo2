from pprint import pprint

from tqdm import tqdm

from chembl_webresource_client.new_client import new_client
from drug_gene import get_tsv_gene_ids
from protein_entrez import get_gene_id

chembl_file = 'interactions.txt'

def get_accession_from_activity(activity):
    return activity[0]['target_components'][0]['accession']


def chembl_id_to_accession_id(chembl_file):
    """
    Returns a dictionary of chembl_id:accession_id pairs, given the name of the interactions file
    @TODO FIND OUT WHY SOME ACCESSION IDS CAN'T BE FOUND
    :param chembl_file: interactions file
    :return: dict
    """
    activities = new_client.target
    accession_ids = dict()
    with open(chembl_file, 'r') as file:
        for line in tqdm(file):
            chembl_id = line.split('\t')[1]
            try:
                res2 = activities.filter(target_chembl_id=[line.strip()],
                                         pchembl_value__isnull=False)
                accession_ids[line.strip()] = get_accession_from_activity(res2)
            except Exception:
                accession_ids[line.strip()] = None
    return accession_ids


accession_dict = chembl_id_to_accession_id(chembl_file)
accession_list = [x for x in list(accession_dict.values()) if x is not None]  # Remove None values

gene_list = get_gene_id(','.join(accession_list))  # Join the accession list into a single string

tsv_genes = get_tsv_gene_ids()  # Get gene ids from tsv

for found_gene in gene_list:  # Compare found genes with tsv genes
    if found_gene in tsv_genes:
        print(found_gene)
