from pprint import pprint

from conda._vendor.tqdm import tqdm

from chembl_webresource_client.new_client import new_client


def get_accession_from_activity(activity):
    return activity[0]['target_components'][0]['accession']


def chembl_id_to_gene_id(chembl_file):
    activities = new_client.target
    accession_ids = dict()
    with open(chembl_file, 'r') as chembl_file_id:
        for line in tqdm(chembl_file_id):
            try:
                res2 = activities.filter(target_chembl_id=[line.strip()],
                                         pchembl_value__isnull=False)
                accession_ids[line.strip()]=get_accession_from_activity(res2)
            except Exception:
                accession_ids[line.strip()] =None
    return accession_ids
