from pprint import pprint

from conda._vendor.tqdm import tqdm

from chembl_webresource_client.new_client import new_client
from drug_gene import get_tsv_gene_ids
from protein_entrez import get_gene_id


def get_accession_from_activity(activity):
    return activity[0]['target_components'][0]['accession']


def chembl_id_to_accession_id(chembl_file):
    """
    Returns a dictionary of chembl_id:accession:id pairs
    :param chembl_file:
    :return:
    """
    activities = new_client.target
    accession_ids = dict()
    with open(chembl_file, 'r') as chembl_file_id:
        for line in tqdm(chembl_file_id):
            try:
                res2 = activities.filter(target_chembl_id=[line.strip()],
                                         pchembl_value__isnull=False)
                accession_ids[line.strip()] = get_accession_from_activity(res2)
            except Exception:
                accession_ids[line.strip()] = None
    return accession_ids


accession_dict = chembl_id_to_accession_id('interactions_sorted.txt')
accession_list = [x for x in list(accession_dict.values()) if x is not None]

gene_list = get_gene_id(','.join(accession_list))

tsv_genes = get_tsv_gene_ids()

for found_gene in gene_list:
    if found_gene in tsv_genes:
        print(found_gene)

